.PHONY: help install dev-install test lint format figures paper-bundle exe clean docker-build docker-run

PYTHON ?= python

help:
	@echo "NLE-Predictor Automation and Research Task Runner"
	@echo "================================================="
	@echo "make install      : Install package in editable mode"
	@echo "make dev-install  : Install package with all development dependencies"
	@echo "make test         : Run complete test suite with pytest"
	@echo "make lint         : Run Ruff code quality and style checks"
	@echo "make format       : Format Python code using Ruff"
	@echo "make figures      : Generate all 600 DPI manuscript publication figures"
	@echo "make paper-bundle : Create self-contained manuscript zip bundle for journal submission"
	@echo "make exe          : Build standalone desktop application executable (.exe)"
	@echo "make clean        : Remove build artifacts, pycache, and temporary files"
	@echo "make docker-build : Build containerized Docker image"
	@echo "make docker-run   : Run containerized CLI"

install:
	$(PYTHON) -m pip install -e .

dev-install:
	$(PYTHON) -m pip install -e ".[dev]"

test:
	$(PYTHON) -m pytest -v

lint:
	$(PYTHON) -m ruff check .

format:
	$(PYTHON) -m ruff format .

figures:
	$(PYTHON) scripts/generate_paper_figures.py

paper-bundle:
	$(PYTHON) -c "import zipfile, os; z = zipfile.ZipFile('manuscript/submission_bundle.zip', 'w', zipfile.ZIP_DEFLATED); [z.write(os.path.join(r, f), os.path.relpath(os.path.join(r, f), 'manuscript')) for r, _, fs in os.walk('manuscript') for f in fs if not f.endswith('.zip')]; print('Created manuscript/submission_bundle.zip')"

exe:
	pyinstaller --clean nle_predictor.spec

clean:
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	rm -rf build/ dist/ .coverage htmlcov/ .mypy_cache/ .ruff_cache/

docker-build:
	docker build -t nle-predictor:latest .

docker-run:
	docker run --rm -it nle-predictor:latest --help

