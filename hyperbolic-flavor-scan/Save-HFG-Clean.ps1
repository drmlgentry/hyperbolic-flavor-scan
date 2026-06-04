# Save-HFG-Clean.ps1
# Run this script from inside C:\dev\hyperbolic-flavor-scan
# It will exclude .venv, __pycache__, .sage-cache, .git, and temporary files.

param(
    [string]$DestinationRoot = "C:\dev\HFG_Backups",
    [switch]$CreateZip,
    [switch]$Verbose
)

$timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
$backupDir = Join-Path $DestinationRoot "HFG_Backup_Clean_$timestamp"
$scriptsDir = Get-Location   # current directory (C:\dev\hyperbolic-flavor-scan)
$papersDir = "C:\dev\framework\papers"
$handoffPath = "C:\dev\framework\HANDOFF_MAY25_2026.md"

Write-Host "Creating clean backup at: $backupDir" -ForegroundColor Cyan
New-Item -ItemType Directory -Path $backupDir -Force | Out-Null

# 1. Copy handoff
if (Test-Path $handoffPath) {
    Copy-Item $handoffPath -Destination $backupDir -Force
    Write-Host "Copied handoff." -ForegroundColor Green
} else {
    Write-Host "Handoff not found." -ForegroundColor Yellow
}

# 2. Copy papers (preserve structure, skip temporary files)
if (Test-Path $papersDir) {
    $papersBackup = Join-Path $backupDir "papers"
    Write-Host "Copying papers from $papersDir ..." -ForegroundColor Cyan
    robocopy $papersDir $papersBackup /E /NJH /NJS /NDL /NFL /XF "*.aux" "*.log" "*.out" "*.toc" "*.synctex.gz" "*.pyc" "*.sage-cache" "*.DS_Store" "__pycache__*" /XD ".git" "__pycache__" ".sage-cache" "venv" ".venv"
    if ($LASTEXITCODE -ge 8) {
        Write-Host "Error during paper copy (code $LASTEXITCODE)" -ForegroundColor Red
    } else {
        Write-Host "Papers copied." -ForegroundColor Green
    }
} else {
    Write-Host "Papers directory not found." -ForegroundColor Yellow
}

# 3. Copy scripts and data (exclude .venv, __pycache__, compiled/temp files)
if (Test-Path $scriptsDir) {
    $scriptsBackup = Join-Path $backupDir "hyperbolic-flavor-scan"
    Write-Host "Copying scripts and data from $scriptsDir ..." -ForegroundColor Cyan
    robocopy $scriptsDir $scriptsBackup /E /NJH /NJS /NDL /NFL /XF "*.aux" "*.log" "*.out" "*.toc" "*.synctex.gz" "*.pyc" "*.sage-cache" "*.DS_Store" "*.json" /XD ".venv" "__pycache__" ".sage-cache" "venv" ".git" "archive" "results" "__pycache__"
    if ($LASTEXITCODE -ge 8) {
        Write-Host "Error during script copy (code $LASTEXITCODE)" -ForegroundColor Red
    } else {
        Write-Host "Scripts and data copied." -ForegroundColor Green
    }
} else {
    Write-Host "Current directory not found." -ForegroundColor Yellow
}

# 4. Optional ZIP
if ($CreateZip) {
    $zipFile = Join-Path $DestinationRoot "HFG_Backup_Clean_$timestamp.zip"
    Write-Host "Creating ZIP archive: $zipFile" -ForegroundColor Cyan
    if (Get-Command Compress-Archive -ErrorAction SilentlyContinue) {
        Compress-Archive -Path $backupDir -DestinationPath $zipFile -Force
        Write-Host "ZIP created." -ForegroundColor Green
    } else {
        Write-Host "Compress-Archive not available; please manually zip." -ForegroundColor Yellow
    }
}

Write-Host "Clean backup completed at $backupDir" -ForegroundColor Cyan