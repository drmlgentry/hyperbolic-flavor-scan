# build-paper.ps1  (cp-admission-prd)
$ROOT   = Split-Path -Parent $PSScriptRoot
$PAPER  = Join-Path $ROOT "paper"
$OUTPUT = Join-Path $ROOT "output"

if (-not (Get-Command pdflatex -ErrorAction SilentlyContinue)) {
    Write-Error "pdflatex not found. TeX Live is installed but not on PATH (or shell not restarted)."
    exit 1
}

# Ensure output directory exists
New-Item -ItemType Directory -Force -Path $OUTPUT | Out-Null

Set-Location $PAPER

# Clean common LaTeX artifacts
Remove-Item *.aux, *.log, *.out, *.toc, *.fls, *.fdb_latexmk -ErrorAction SilentlyContinue

Write-Host "Compiling (pass 1)..."
pdflatex -interaction=nonstopmode main.tex

Write-Host "Compiling (pass 2)..."
pdflatex -interaction=nonstopmode main.tex

# Move PDF to output
Move-Item -Force main.pdf (Join-Path $OUTPUT "main.pdf")

Write-Host "Done: $OUTPUT\main.pdf"
