# syntax=docker/dockerfile:1
FROM python:3.11-slim as base

# Set environment variables for reproducibility and headless execution
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    PIP_NO_CACHE_DIR=1

# Install essential system libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libgomp1 \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy dependency specifications first to leverage Docker layer caching
COPY pyproject.toml requirements.txt requirements-dev.txt ./

# Install dependencies and package in editable mode
COPY src/ ./src/
COPY inputs/ ./inputs/
COPY examples/ ./examples/
COPY tests/ ./tests/
COPY README.md LICENSE ./

RUN pip install --upgrade pip setuptools wheel \
    && pip install -e ".[dev]"

# Ensure unprivileged user execution for security and academic HPC clusters
RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /app
USER appuser

WORKDIR /workspace

# Default entrypoint runs the unified CLI
ENTRYPOINT ["nle-predictor"]
CMD ["--help"]
