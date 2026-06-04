# Save-HFGState.ps1
# Run this in PowerShell to archive all HFG project files.
# It will create a timestamped folder (e.g., HFG_Backup_20260526_143022) and copy:
#   - All paper LaTeX files and PDFs from C:\dev\framework\papers
#   - All scripts and data from C:\dev\hyperbolic-flavor-scan
#   - The handoff markdown file
# Optionally, it will also create a ZIP archive.

param(
    [string]$DestinationRoot = "C:\dev\HFG_Backups",
    [switch]$CreateZip
)

$timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
$backupDir = Join-Path $DestinationRoot "HFG_Backup_$timestamp"
$scriptsDir = "C:\dev\hyperbolic-flavor-scan"
$papersDir = "C:\dev\framework\papers"
$handoffPath = "C:\dev\framework\HANDOFF_MAY25_2026.md"

Write-Host "Creating backup directory: $backupDir" -ForegroundColor Cyan
New-Item -ItemType Directory -Path $backupDir -Force | Out-Null

# Copy handoff
if (Test-Path $handoffPath) {
    Copy-Item $handoffPath -Destination $backupDir -Force
    Write-Host "Copied handoff document." -ForegroundColor Green
} else {
    Write-Host "Handoff document not found at $handoffPath" -ForegroundColor Yellow
}

# Copy papers
if (Test-Path $papersDir) {
    $papersBackup = Join-Path $backupDir "papers"
    New-Item -ItemType Directory -Path $papersBackup -Force | Out-Null
    # Recursively copy all .tex, .pdf, .bib, .cls, .sty, .fig, .png, .jpg
    Get-ChildItem -Path $papersDir -Recurse -Include *.tex,*.pdf,*.bib,*.cls,*.sty,*.fig,*.png,*.jpg,*.md |
        Copy-Item -Destination { Join-Path $papersBackup $_.FullName.Substring($papersDir.Length) } -Force
    Write-Host "Copied papers directory." -ForegroundColor Green
} else {
    Write-Host "Papers directory not found at $papersDir" -ForegroundColor Yellow
}

# Copy scripts and data (excluding large temporary or compiled files)
if (Test-Path $scriptsDir) {
    $scriptsBackup = Join-Path $backupDir "hyperbolic-flavor-scan"
    New-Item -ItemType Directory -Path $scriptsBackup -Force | Out-Null
    # Exclude common temp files, .pyc, __pycache__, .sage-cache, .ipynb_checkpoints, and large data dumps
    Get-ChildItem -Path $scriptsDir -Recurse |
        Where-Object {
            $_.FullName -notmatch "__pycache__" -and
            $_.FullName -notmatch "\.pyc$" -and
            $_.FullName -notmatch "\.sage-cache" -and
            $_.FullName -notmatch "\.ipynb_checkpoints" -and
            $_.FullName -notmatch "\.DS_Store" -and
            $_.FullName -notmatch "\.log$" -and
            $_.Extension -notin ".aux",".log",".out",".toc",".synctex.gz"
        } |
        Copy-Item -Destination { Join-Path $scriptsBackup $_.FullName.Substring($scriptsDir.Length) } -Force
    Write-Host "Copied scripts and data directory." -ForegroundColor Green
} else {
    Write-Host "Scripts directory not found at $scriptsDir" -ForegroundColor Yellow
}

# Optionally create a ZIP archive of the backup folder
if ($CreateZip) {
    $zipFile = Join-Path $DestinationRoot "HFG_Backup_$timestamp.zip"
    Write-Host "Creating ZIP archive: $zipFile" -ForegroundColor Cyan
    if (Get-Command Compress-Archive -ErrorAction SilentlyContinue) {
        Compress-Archive -Path $backupDir -DestinationPath $zipFile -Force
        Write-Host "ZIP archive created." -ForegroundColor Green
    } else {
        Write-Host "Compress-Archive not available; please manually ZIP the folder." -ForegroundColor Yellow
    }
}

Write-Host "Backup completed at $backupDir" -ForegroundColor Cyan