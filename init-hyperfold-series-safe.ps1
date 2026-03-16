Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# ---- CONFIG ----
$Root = "C:\dev\hyperfold-series"

# ---- SAFETY ----
if ([string]::IsNullOrWhiteSpace($Root)) { throw "Root is empty. Aborting." }
$Root = [System.IO.Path]::GetFullPath($Root)

$forbidden = @("C:\", "C:\dev")
if ($forbidden -contains $Root.TrimEnd("\")) { throw "Refusing forbidden Root=$Root" }
if (-not ($Root.ToLower().StartsWith("c:\dev\"))) { throw "Refusing Root outside C:\dev: $Root" }

# ---- CREATE DIRS ----
$dirs = @(
  $Root,
  (Join-Path $Root "docs"),
  (Join-Path $Root "spine"),
  (Join-Path $Root "papers"),
  (Join-Path $Root "papers\paper1_geometry"),
  (Join-Path $Root "papers\paper2_certification"),
  (Join-Path $Root "papers\paper3_phenomenology"),
  (Join-Path $Root "papers\paper4_overlap_framework"),
  (Join-Path $Root "papers\paper5_dark_sector_optional"),
  (Join-Path $Root "scripts"),
  (Join-Path $Root "imports")
)
foreach ($d in $dirs) { New-Item -ItemType Directory -Force -Path $d | Out-Null }

# ---- FILES ----
@"
# LaTeX build
*.aux
*.bbl
*.blg
*.log
*.out
*.toc
*.lof
*.lot
*.fls
*.fdb_latexmk
*.synctex.gz
build/
results/
"@ | Set-Content -Encoding UTF8 (Join-Path $Root ".gitignore")

@"
* text=auto
*.tex text eol=lf
*.bib text eol=lf
*.py  text eol=lf
*.ps1 text eol=crlf
"@ | Set-Content -Encoding UTF8 (Join-Path $Root ".gitattributes")

@"
# hyperfold-series

Private hub repo for the consolidated program.

- `spine/` : dependency skeleton (canonical)
- `papers/`: composed papers with provenance
- `imports/`: provenance logs + source hashes
"@ | Set-Content -Encoding UTF8 (Join-Path $Root "README.md")

@"
# Dependency Skeleton (Program Spine)
See spine/spine.tex for the compilable LaTeX version.
"@ | Set-Content -Encoding UTF8 (Join-Path $Root "docs\DEPENDENCY_SKELETON.md")

@"
\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb,amsthm}
\title{Program Spine: Dependency Skeleton}
\author{Marvin Gentry}
\date{\today}
\begin{document}\maketitle
% Full skeleton text will be inserted here.
\end{document}
"@ | Set-Content -Encoding UTF8 (Join-Path $Root "spine\spine.tex")

@"
Set-StrictMode -Version Latest
`$ErrorActionPreference = "Stop"
Push-Location (Split-Path -Parent `$MyInvocation.MyCommand.Path)
try {
  Push-Location "spine"
  if (-not (Get-Command latexmk -ErrorAction SilentlyContinue)) {
    throw "latexmk not found. Install TeX Live/MiKTeX with latexmk."
  }
  latexmk -C | Out-Null
  latexmk -pdf -interaction=nonstopmode -halt-on-error spine.tex
  Write-Host "OK: spine\spine.pdf"
} finally {
  Pop-Location
  Pop-Location
}
"@ | Set-Content -Encoding UTF8 (Join-Path $Root "build.ps1")

# ---- INIT GIT (ONLY INSIDE ROOT) ----
Push-Location $Root
try {
  if (-not (Test-Path -LiteralPath ".\.git")) {
    git init | Out-Null
    git add . | Out-Null
    git commit -m "Initialize hyperfold-series hub repo" | Out-Null
  }
  Write-Host "OK: initialized at $Root"
} finally {
  Pop-Location
}

Write-Host "DONE. Repo created at: $Root"