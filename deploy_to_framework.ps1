# =============================================================================
# deploy_to_framework.ps1
# Deploys files from Downloads folder to their correct locations in C:\dev\framework
# Usage: Run from PowerShell. Adjust $DownloadsDir if needed.
# =============================================================================

$DownloadsDir = "$env:USERPROFILE\Downloads"
$FrameworkDir = "C:\dev\framework"
$ScanDir      = "C:\dev\hyperbolic-flavor-scan"

# -----------------------------------------------------------------------------
# FILE ROUTING TABLE
# Add new files here as needed: @{File="filename.ext"; Dest="relative\path"}
# -----------------------------------------------------------------------------
$routes = @(
    # AHP submission -- corrected CP holonomy paper
    @{ File = "gentry-holonomy-cp-ahp.tex";
       Dest = "$FrameworkDir\papers\holonomy-cp" },

    # Lucas structure paper (reframed)
    @{ File = "gentry_lucas_structure.tex";
       Dest = "$FrameworkDir\papers\lucas-structure" },

    # Weeks-Dehn paper (CS corrected)
    @{ File = "gentry-weeks-dehn.tex";
       Dest = "$FrameworkDir\papers\weeks-dehn" },

    # Hyperbolic lattice paper (CS corrected)
    @{ File = "gentry-hyperbolic-lattice.tex";
       Dest = "$FrameworkDir\papers\hyperbolic-lattice" },

    # Reproducibility script
    @{ File = "hfg_reproduce.py";
       Dest = $ScanDir }
)

# -----------------------------------------------------------------------------
Write-Host "=== HFG Framework Deployment Script ===" -ForegroundColor Cyan
Write-Host "Downloads: $DownloadsDir"
Write-Host "Framework: $FrameworkDir"
Write-Host ""

$deployed = 0
$skipped  = 0
$missing  = 0

foreach ($route in $routes) {
    $src  = Join-Path $DownloadsDir $route.File
    $dest = $route.Dest

    if (-not (Test-Path $src)) {
        Write-Host "  MISSING : $($route.File)" -ForegroundColor Yellow
        $missing++
        continue
    }

    # Create destination directory if needed
    if (-not (Test-Path $dest)) {
        New-Item -ItemType Directory -Path $dest -Force | Out-Null
        Write-Host "  CREATED : $dest" -ForegroundColor DarkGray
    }

    $destFile = Join-Path $dest $route.File
    Copy-Item -Path $src -Destination $destFile -Force
    Write-Host "  DEPLOYED: $($route.File) -> $dest" -ForegroundColor Green
    $deployed++
}

Write-Host ""
Write-Host "Done: $deployed deployed, $skipped skipped, $missing missing." -ForegroundColor Cyan

# -----------------------------------------------------------------------------
# Git commit deployed files
# -----------------------------------------------------------------------------
Write-Host ""
$doCommit = Read-Host "Commit deployed files to git? (y/n)"
if ($doCommit -eq 'y') {
    Set-Location $FrameworkDir
    git add papers\holonomy-cp\gentry-holonomy-cp-ahp.tex `
            papers\lucas-structure\gentry_lucas_structure.tex `
            papers\weeks-dehn\gentry-weeks-dehn.tex `
            papers\hyperbolic-lattice\gentry-hyperbolic-lattice.tex 2>$null
    git commit -m "Deploy corrected files: AHP CP holonomy (email+funding+RINP refs), Lucas reframe, CS corrections"
    git push origin main
    Write-Host "Committed and pushed." -ForegroundColor Green
}
