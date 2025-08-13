#!/bin/bash
# Setup GitHub repository

echo "Initializing Git repository..."
git init
git add .
git commit -m "Initial commit: Code and documentation for GBM UPR Rheostat analysis"

echo "Repository ready!"
echo "Next steps:"
echo "1. Create repository on GitHub: https://github.com/new"
echo "2. Add remote: git remote add origin https://github.com/zduzgun/GBM-UPR-Rheostat.git"
echo "3. Push: git push -u origin main"
