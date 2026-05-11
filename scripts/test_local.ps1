# Run before every push to verify code quality and tests
# Usage: .\scripts\test_local.ps1 [-SkipFormat]
param([switch]$SkipFormat)

$ErrorActionPreference = "Stop"

if (-not $SkipFormat) {
    Write-Host "=== Runic formatting check ===" -ForegroundColor Cyan
    julia --project=@runic -e 'using Runic; exit(Runic.main(["--check", "src/", "test/", "ext/"]))'
    if ($LASTEXITCODE -ne 0) {
        Write-Host "Formatting issues found. Auto-fixing..." -ForegroundColor Yellow
        julia --project=@runic -e 'using Runic; Runic.main(["--inplace", "src/", "test/", "ext/"])'
        Write-Host "Fixed. Please review and stage the changes." -ForegroundColor Yellow
        exit 1
    }
}

Write-Host "=== Running tests ===" -ForegroundColor Cyan
julia --project -e 'using Pkg; Pkg.test()'
if ($LASTEXITCODE -ne 0) {
    Write-Host "Tests FAILED." -ForegroundColor Red
    exit 1
}

Write-Host "All checks passed. Ready to push." -ForegroundColor Green
