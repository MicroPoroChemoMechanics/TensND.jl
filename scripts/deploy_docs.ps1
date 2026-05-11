# Build and deploy documentation to Codeberg Pages (gh-pages branch)
# Prerequisites: DOCUMENTER_KEY must be set (see WORKFLOW.md for one-time setup)
# Usage: .\scripts\deploy_docs.ps1

$ErrorActionPreference = "Stop"

if (-not $env:DOCUMENTER_KEY) {
    Write-Error @"
DOCUMENTER_KEY is not set.
Run the one-time setup from WORKFLOW.md (etape 1 a 3) :
  [PowerShell] ssh-keygen -t ed25519 -C "documenter@codeberg" -f "$env:USERPROFILE\.ssh\documenter_codeberg" -N ''
  [PowerShell] `$bytes = [System.IO.File]::ReadAllBytes("$env:USERPROFILE\.ssh\documenter_codeberg")
  [PowerShell] `$b64 = [Convert]::ToBase64String(`$bytes)
  [PowerShell] [System.Environment]::SetEnvironmentVariable("DOCUMENTER_KEY", `$b64, "User")
Puis ouvre un nouveau terminal et relance ce script.
"@
    exit 1
}

Write-Host "=== Running doctests ===" -ForegroundColor Cyan
julia --project=docs/ -e '
    using Documenter: DocMeta, doctest
    using TensND
    DocMeta.setdocmeta!(TensND, :DocTestSetup,
        :(using TensND, LinearAlgebra, SymPy, Tensors, OMEinsum, Rotations);
        recursive=true)
    doctest(TensND)
'
if ($LASTEXITCODE -ne 0) {
    Write-Host "Doctests FAILED." -ForegroundColor Red
    exit 1
}

Write-Host "=== Building and deploying documentation ===" -ForegroundColor Cyan
$env:CI = "true"
julia --project=docs/ docs/make.jl
$env:CI = $null

if ($LASTEXITCODE -ne 0) {
    Write-Host "Documentation build FAILED." -ForegroundColor Red
    exit 1
}

Write-Host "Documentation deployed to Codeberg Pages." -ForegroundColor Green
Write-Host "Dev:    https://MicroPoroChemoMechanics.codeberg.page/TensND.jl/dev/" -ForegroundColor Cyan
Write-Host "Stable: https://MicroPoroChemoMechanics.codeberg.page/TensND.jl/stable/" -ForegroundColor Cyan
