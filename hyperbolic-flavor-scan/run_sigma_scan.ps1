$sigmas = @(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
$working_script = "scan_flavors_working.py"
$jobs = @()

foreach ($s in $sigmas) {
    $outfile = "scan_results_sigma_$s.csv"
    $job = Start-Job -ScriptBlock {
        param($s, $outfile, $script)
        Set-Location $using:PWD
        python $script $s $outfile
    } -ArgumentList $s, $outfile, $working_script
    $jobs += $job
}

Write-Host "Waiting for all sigma scans to complete..."
$jobs | Wait-Job | Receive-Job
Write-Host "All done!"