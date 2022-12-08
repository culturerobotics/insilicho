install:
	pip install --upgrade pip
	pip install poetry==1.1.12
ifeq ($(shell uname -sm),Darwin arm64)
	OPENBLAS="$(shell brew --prefix openblas)" poetry install --no-root ${POETRY_ARGS}
else
	poetry install --no-root ${POETRY_ARGS}
endif

tc typecheck:
	poetry run mypy . --show-error-codes

test: tc
	poetry run pytest --timeout=5 --cov-report=xml --cov-report=html --cov=. --junitxml=test-metadata/junit.xml

test-verbose: tc
	poetry run pytest -vvv --capture=no --timeout=5 --cov-report=xml --cov-report=html --cov=. --junitxml=test-metadata/junit.xml

.PHONY: docs-local

docs-local:
	cd docs; make clean && make html; open _build/html/index.html