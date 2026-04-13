install:
	uv sync
	uv run pre-commit install

tc typecheck:
	uv run mypy . --show-error-codes

test: tc
	uv run pytest --timeout=5 --cov-report=xml --cov-report=html --cov=. --junitxml=test-metadata/junit.xml -vvv

test-verbose: tc
	uv run pytest -vvv --capture=no --timeout=5 --cov-report=xml --cov-report=html --cov=. --junitxml=test-metadata/junit.xml

.PHONY: docs-local

docs-local:
	cd docs; SPHINXBUILD="uv run sphinx-build" make clean && SPHINXBUILD="uv run sphinx-build" make html; open _build/html/index.html
