# Save-HFGState-Robust.ps1
param(
    [string]$DestinationRoot = "C:\dev\HFG_Backups",
    [switch]$CreateZip,
    [switch]$Verbose
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
    Write-Host "Handoff document not found." -ForegroundColor Yellow
}

# Function to copy directory with structure, skipping problematic files
function Copy-DirectoryStructure {
    param([string]$Source, [string]$Dest)
    if (-not (Test-Path $Source)) {
        Write-Host "Source directory not found: $Source" -ForegroundColor Yellow
        return
    }
    Get-ChildItem -Path $Source -Recurse -ErrorAction SilentlyContinue | ForEach-Object {
        # Skip files with invalid characters or names (e.g., "1$.")
        if ($_.FullName -match '[\$]' -or $_.Name -match '^\d+\$\.') {
            Write-Host "  Skipping invalid file: $($_.FullName)" -ForegroundColor DarkYellow
            return
        }
        $relativePath = $_.FullName.Substring($Source.Length + 1)
        $destPath = Join-Path $Dest $relativePath
        if ($_.PSIsContainer) {
            New-Item -ItemType Directory -Path $destPath -Force -ErrorAction SilentlyContinue | Out-Null
        } else {
            # Filter unwanted extensions
            $ext = $_.Extension.ToLower()
            if ($ext -in '.aux','.log','.out','.toc','.synctex.gz','.pyc','.sage-cache','.DS_Store') { return }
            if ($_.Name -like '__pycache__*') { return }
            try {
                Copy-Item -Path $_.FullName -Destination $destPath -Force -ErrorAction Stop
                if ($Verbose) { Write-Host "  Copied $relativePath" }
            } catch {
                Write-Host "  Failed to copy $relativePath : $_" -ForegroundColor Red
            }
        }
    }
}

# Copy papers
$papersBackup = Join-Path $backupDir "papers"
New-Item -ItemType Directory -Path $papersBackup -Force | Out-Null
Write-Host "Copying papers from $papersDir ..." -ForegroundColor Cyan
Copy-DirectoryStructure -Source $papersDir -Dest $papersBackup
Write-Host "Finished copying papers." -ForegroundColor Green

# Copy scripts and data
$scriptsBackup = Join-Path $backupDir "hyperbolic-flavor-scan"
New-Item -ItemType Directory -Path $scriptsBackup -Force | Out-Null
Write-Host "Copying scripts and data from $scriptsDir ..." -ForegroundColor Cyan
Copy-DirectoryStructure -Source $scriptsDir -Dest $scriptsBackup
Write-Host "Finished copying scripts." -ForegroundColor Green

# Optionally create ZIP
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