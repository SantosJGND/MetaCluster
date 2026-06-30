# Contributing

Thank you for your interest in contributing to the Metagenomics Evaluation Pipeline.

## Development Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/insa/metagenomics-evaluation-pipeline.git
   cd metagenomics-evaluation-pipeline
   ```

2. Create and activate a virtual environment:
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   ```

3. Install in development mode:
   ```bash
   pip install -e .
   ```

4. Install pre-commit hooks:
   ```bash
   pre-commit install
   ```

## Running Tests

Run all tests:
```bash
pytest
```

Run unit tests only (skip integration, slow, and NCBI-dependent tests):
```bash
pytest -m "not slow and not integration and not requires_ncbi"
```

Run with coverage:
```bash
pytest --cov
```

## Code Style

This project uses [ruff](https://docs.astral.sh/ruff/) for linting and formatting:

```bash
ruff check .     # lint
ruff format .    # format
```

Pre-commit hooks will enforce style automatically.

## Pull Request Process

1. Create a feature branch from `main`
2. Make your changes
3. Ensure all tests pass
4. Run the linter and formatter
5. Submit a PR with a clear description of the changes
