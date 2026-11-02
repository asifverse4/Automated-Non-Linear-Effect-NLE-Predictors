# PowerShell Task Runner for NLE-Predictor (Windows)
param (
    [Parameter(Position=0)]
    [string]$Target = "help"
)

function Show-Help {
    Write-Host "NLE-Predictor PowerShell Task Runner" -ForegroundColor Cyan
    Write-Host "====================================" -ForegroundColor Cyan
    Write-Host ".\make.ps1 install     : Install package in editable mode"
    Write-Host ".\make.ps1 dev-install : Install package with dev dependencies"
    Write-Host ".\make.ps1 test        : Run complete test suite with pytest"
    Write-Host ".\make.ps1 lint        : Run Ruff linter"
    Write-Host ".\make.ps1 format      : Format code with Ruff"
    Write-Host ".\make.ps1 figures      : Regenerate manuscript figures"
    Write-Host ".\make.ps1 paper-bundle : Create self-contained manuscript zip bundle for submission"
    Write-Host ".\make.ps1 exe          : Build standalone Windows desktop application (.exe)"
    Write-Host ".\make.ps1 clean        : Remove caches and build artifacts"
}

switch ($Target.ToLower()) {
    "help" { Show-Help }
    "install" { python -m pip install -e . }
    "dev-install" { python -m pip install -e ".[dev]" }
    "test" { python -m pytest -v }
    "lint" { python -m ruff check . }
    "format" { python -m ruff format . }
    "figures" { python scripts/generate_paper_figures.py }
    "paper-bundle" {
        python -c "import zipfile, os; z = zipfile.ZipFile('manuscript/submission_bundle.zip', 'w', zipfile.ZIP_DEFLATED); [z.write(os.path.join(r, f), os.path.relpath(os.path.join(r, f), 'manuscript')) for r, _, fs in os.walk('manuscript') for f in fs if not f.endswith('.zip')]; print('Created manuscript/submission_bundle.zip')"
    }
    "exe" {
        pyinstaller --clean nle_predictor.spec
        Compress-Archive -Path dist/nle-predictor -DestinationPath dist/nle-predictor-windows-x64.zip -Force
        Write-Host "Standalone executable packaged: dist/nle-predictor-windows-x64.zip" -ForegroundColor Green
    }
    "clean" {
        Get-ChildItem -Recurse -Directory -Filter "__pycache__" | Remove-Item -Recurse -Force -ErrorAction SilentlyContinue
        Get-ChildItem -Recurse -Directory -Filter ".pytest_cache" | Remove-Item -Recurse -Force -ErrorAction SilentlyContinue
        Get-ChildItem -Recurse -Directory -Filter "*.egg-info" | Remove-Item -Recurse -Force -ErrorAction SilentlyContinue
        Remove-Item -Recurse -Force build, dist, .mypy_cache, .ruff_cache -ErrorAction SilentlyContinue
        Write-Host "Clean completed." -ForegroundColor Green
    }
    default {
        Write-Host "Unknown target: $Target" -ForegroundColor Red
        Show-Help
    }
}

