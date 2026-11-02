.PHONY: help install dev-install test lint format figures clean docker-build docker-run

PYTHON ?= python

help:
	@echo "NLE-Predictor Automation and Research Task Runner"
	@echo "================================================="
	@echo "make install      : Install package in editable mode"
	@echo "make dev-install  : Install package with all development dependencies"
	@echo "make test         : Run complete 49-test suite with pytest"
	@echo "make lint         : Run Ruff code quality and style checks"
	@echo "make format       : Format Python code using Ruff"
	@echo "make figures      : Generate all 600 DPI manuscript publication figures"
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

clean:
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	rm -rf build/ dist/ .coverage htmlcov/ .mypy_cache/ .ruff_cache/

docker-build:
	docker build -t nle-predictor:latest .

docker-run:
	docker run --rm -it nle-predictor:latest --help
