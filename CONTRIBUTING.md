# Contributing

## Development setup
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .[dev]
pre-commit install
```

## Before opening a PR
Run all local checks:
```bash
ruff check .
python -m unittest discover -s tests -p "test_*.py" -v
bandit -r . -x tests,.venv
pip-audit
```

## Pull request requirements
- Keep changes scoped and include tests for behavior changes.
- Do not commit secrets or credentials.
- Ensure CI and Security workflows are green.
