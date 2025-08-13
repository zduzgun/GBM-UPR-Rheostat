# ===================================================================
# KALLISTO RESULTS AGGREGATION - ROBUST PARSER
# ===================================================================

# --- Load Packages ---
suppressMessages({
  library(tximport)
  library(readr)
  library(dplyr)
})

options(readr.num_columns = 0, stringsAsFactors = FALSE)
gc()

cat("?? ROBUST PARSER: KALLISTO AGGREGATION\n")
cat("======================================\n")

# --- CONFIGURATION ---
CONFIG <- list(
  kallisto_dir = "kallisto_output",
  gencode_fasta = "referans/gencode.v46.transcripts.fa.gz",
  output_file = "glioblastoma_gene_counts.csv",
  chunk_size = 30,
  t2g_cache_file = "t2g_mapping_robust.rds"
)

# --- ROBUST FASTA PARSING ---
create_robust_gencode_parser <- function(fasta_path, cache_file = NULL) {
  
  # Check cache
  if(!is.null(cache_file) && file.exists(cache_file)) {
    cat("?? Loading cached t2g mapping...\n")
    cached <- readRDS(cache_file)
    cat(paste("? Loaded", nrow(cached), "cached mappings\n"))
    return(cached)
  }

  cat("?? Parsing GENCODE FASTA with robust method...\n")
  
  con <- file(fasta_path, "r")
  
  # Test different splitting methods with first header
  test_line <- readLines(con, n = 1000)
  first_header <- test_line[startsWith(test_line, ">")][1]
  
  cat("?? Testing splitting methods:\n")
  cat(paste("Sample header:", substr(first_header, 1, 100), "...\n"))
  
  clean_header <- gsub("^>", "", first_header)
  
  # Try different splitting approaches
  splitting_methods <- list(
    # Method 1: Standard pipe
    method1 = function(x) strsplit(x, "\\|", fixed = TRUE)[[1]],
    
    # Method 2: Different pipe patterns
    method2 = function(x) strsplit(x, "|", fixed = TRUE)[[1]],
    
    # Method 3: Unicode pipe variations
    method3 = function(x) strsplit(x, "\u007C")[[1]],  # Unicode pipe
    
    # Method 4: Regex approach
    method4 = function(x) unlist(strsplit(x, "[||]")),  # Regular and fullwidth pipe
    
    # Method 5: Manual character splitting
    method5 = function(x) {
      # Find all pipe positions manually
      chars <- strsplit(x, "")[[1]]
      pipe_positions <- which(chars %in% c("|", "|", "\u007C"))
      
      if(length(pipe_positions) == 0) return(x)
      
      # Split manually
      parts <- c()
      start_pos <- 1
      
      for(pos in pipe_positions) {
        parts <- c(parts, substr(x, start_pos, pos - 1))
        start_pos <- pos + 1
      }
      # Add last part
      parts <- c(parts, substr(x, start_pos, nchar(x)))
      
      return(parts)
    }
  )
  
  working_method <- NULL
  working_method_name <- NULL
  
  for(method_name in names(splitting_methods)) {
    cat(paste("   Testing", method_name, "... "))
    
    tryCatch({
      parts <- splitting_methods[[method_name]](clean_header)
      if(length(parts) > 5) {
        cat(paste("? SUCCESS -", length(parts), "parts\n"))
        cat(paste("      Parts: [1]", parts[1], "[6]", parts[6], "\n"))
        working_method <- splitting_methods[[method_name]]
        working_method_name <- method_name
        break
      } else {
        cat(paste("?", length(parts), "parts\n"))
      }
    }, error = function(e) {
      cat("? Error\n")
    })
  }
  
  if(is.null(working_method)) {
    stop("? No splitting method worked! Check FASTA format manually.")
  }
  
  cat(paste("?? Using", working_method_name, "for parsing\n\n"))
  
  # Reset connection
  close(con)
  con <- file(fasta_path, "r")
  
  # Parse entire file with working method
  all_mappings <- list()
  chunk_counter <- 1
  total_processed <- 0
  
  while(TRUE) {
    lines <- readLines(con, n = 50000)
    if(length(lines) == 0) break
    
    headers <- lines[startsWith(lines, ">")]
    if(length(headers) == 0) next
    
    total_processed <- total_processed + length(headers)
    cat(paste("\r   Processing chunk", chunk_counter, "- Total headers:", total_processed))
    
    for(header in headers) {
      clean_header <- gsub("^>", "", header)
      
      tryCatch({
        parts <- working_method(clean_header)
        
        if(length(parts) >= 6) {
          transcript_id <- parts[1]
          gene_id <- if(length(parts) >= 2) parts[2] else NA
          gene_name <- parts[6]
          gene_type <- if(length(parts) >= 8) parts[8] else NA
          
          # Basic quality check
          if(!is.na(gene_name) && gene_name != "" && gene_name != "-" && 
             nchar(gene_name) > 0) {
            
            all_mappings[[length(all_mappings) + 1]] <- list(
              transcript_id = transcript_id,
              gene_id = gene_id, 
              gene_name = gene_name,
              gene_type = gene_type
            )
          }
        }
      }, error = function(e) {
        # Skip problematic headers
      })
    }
    
    chunk_counter <- chunk_counter + 1
    if(chunk_counter %% 20 == 0) gc()
  }
  
  close(con)
  cat("\n")
  
  if(length(all_mappings) == 0) {
    stop("? No valid mappings extracted!")
  }
  
  # Convert to dataframe
  t2g_df <- data.frame(
    transcript_id = sapply(all_mappings, function(x) x$transcript_id),
    gene_id = sapply(all_mappings, function(x) x$gene_id),
    gene_name = sapply(all_mappings, function(x) x$gene_name),
    gene_type = sapply(all_mappings, function(x) x$gene_type),
    stringsAsFactors = FALSE
  )
  
  cat(paste("?? Raw mappings extracted:", nrow(t2g_df), "\n"))
  
  # Apply smart filtering
  cat("?? Applying intelligent filtering...\n")
  
  # Analyze gene_name patterns first
  gene_name_types <- table(t2g_df$gene_name)
  common_types <- names(gene_name_types)[gene_name_types > 1000]
  
  cat("Most common 'gene names' (likely types to filter):\n")
  print(head(sort(gene_name_types, decreasing = TRUE), 10))
  
  # Smart filtering - remove obvious types but keep real genes
  filtered_t2g <- t2g_df %>%
    filter(
      # Basic cleanup
      !is.na(gene_name),
      gene_name != "",
      gene_name != "-",
      nchar(gene_name) > 0,
      
      # Remove common RNA types
      !gene_name %in% c("lncRNA", "miRNA", "snoRNA", "snRNA", "rRNA", 
                       "scRNA", "Mt_tRNA", "ribozyme"),
      
      # Remove obvious pseudogene markers
      !grepl("pseudogene", gene_name, ignore.case = TRUE),
      
      # Remove Ensembl IDs that leaked into gene_name
      !startsWith(gene_name, "ENSG"),
      !startsWith(gene_name, "ENST"),
      
      # Remove very generic names
      !gene_name %in% c("protein_coding", "processed_transcript", 
                       "antisense", "sense_intronic", "TEC")
    )
  
  cat(paste("?? After filtering:", nrow(filtered_t2g), "high-quality mappings\n"))
  
  # Remove duplicates
  final_t2g <- filtered_t2g %>%
    distinct(transcript_id, .keep_all = TRUE)
  
  cat(paste("? Final unique mappings:", nrow(final_t2g), "\n"))
  
  # Show examples
  cat("?? Sample gene mappings:\n")
  print(head(final_t2g %>% select(transcript_id, gene_name), 10))
  
  # Cache result
  if(!is.null(cache_file)) {
    cat("?? Saving to cache...\n")
    saveRDS(final_t2g, cache_file)
  }
  
  return(final_t2g)
}

# --- MAIN EXECUTION ---

# 1. Parse GENCODE with robust method
t2g_mapping <- create_robust_gencode_parser(CONFIG$gencode_fasta, CONFIG$t2g_cache_file)

# 2. Find Kallisto files
cat("?? Finding Kallisto files...\n")
files <- list.files(
  path = CONFIG$kallisto_dir,
  pattern = "abundance\\.h5$", 
  recursive = TRUE,
  full.names = TRUE
)

if(length(files) == 0) {
  stop(paste("? No Kallisto files found in", CONFIG$kallisto_dir))
}

sample_names <- basename(dirname(files))
names(files) <- sample_names

cat(paste("? Found", length(files), "samples\n"))

# 3. Test format compatibility
cat("?? Testing t2g compatibility with Kallisto files...\n")

test_formats <- list(
  # Format 1: Direct match
  direct = t2g_mapping %>% select(transcript_id, gene_name),
  
  # Format 2: Remove version numbers  
  no_version = t2g_mapping %>%
    mutate(transcript_id = gsub("\\.\\d+$", "", transcript_id)) %>%
    select(transcript_id, gene_name) %>%
    distinct(transcript_id, .keep_all = TRUE)
)

working_format <- NULL

for(format_name in names(test_formats)) {
  cat(paste("   Testing", format_name, "format ... "))
  
  tryCatch({
    test_result <- tximport(
      files[1],
      type = "kallisto",
      tx2gene = test_formats[[format_name]],
      ignoreTxVersion = TRUE,
      ignoreAfterBar = TRUE
    )
    
    genes_found <- nrow(test_result$counts)
    if(genes_found > 100) {
      cat(paste("? SUCCESS -", genes_found, "genes\n"))
      working_format <- list(
        t2g = test_formats[[format_name]],
        name = format_name
      )
      rm(test_result)
      break
    }
    rm(test_result)
    cat("? Too few genes\n")
    
  }, error = function(e) {
    cat("? Error\n")
  })
}

if(is.null(working_format)) {
  stop("? No compatible format found!")
}

cat(paste("?? Using format:", working_format$name, "\n\n"))

# 4. Process all samples
total_samples <- length(files)
file_chunks <- split(files, ceiling(seq_along(files) / CONFIG$chunk_size))

cat(paste("? Processing", total_samples, "samples in", length(file_chunks), "chunks\n"))

all_results <- list()
successful <- 0

for(i in seq_along(file_chunks)) {
  cat(paste("?? Chunk", i, "/", length(file_chunks), "... "))
  
  tryCatch({
    txi <- tximport(
      file_chunks[[i]],
      type = "kallisto", 
      tx2gene = working_format$t2g,
      ignoreTxVersion = TRUE,
      ignoreAfterBar = TRUE
    )
    
    all_results[[i]] <- as.data.frame(txi$counts)
    successful <- successful + 1
    cat(paste("?", nrow(txi$counts), "genes\n"))
    rm(txi)
    
  }, error = function(e) {
    cat("? Failed\n")
  })
}

# 5. Combine and save
cat(paste("\n?? Successfully processed:", successful, "/", length(file_chunks), "chunks\n"))

if(successful == 0) {
  stop("? No chunks succeeded!")
}

# Combine results
final_matrix <- all_results[[which(!sapply(all_results, is.null))[1]]]
successful_indices <- which(!sapply(all_results, is.null))

if(length(successful_indices) > 1) {
  for(i in successful_indices[-1]) {
    final_matrix <- cbind(final_matrix, all_results[[i]])
  }
}

# Save
cat("?? Saving results...\n")
write.csv(final_matrix, CONFIG$output_file, quote = FALSE, row.names = TRUE)

cat("\n?? SUCCESS! ??\n")
cat(paste("?? Output:", CONFIG$output_file, "\n"))  
cat(paste("?? Genes:", nrow(final_matrix), "\n"))
cat(paste("?? Samples:", ncol(final_matrix), "\n"))