.PHONY: help install test lint sweep figures paper clean

help:
	@echo "install   editable install with dev extras"
	@echo "test      run the test suite"
	@echo "lint      run ruff and the prose check"
	@echo "sweep     run the full detectability sweep"
	@echo "figures   draw the figures from results/sweep.csv"
	@echo "paper     build paper/main.pdf and the DOCX"
	@echo "clean     remove build artefacts"

install:
	python -m pip install -e ".[dev]"

test:
	pytest -q

lint:
	ruff check src tests
	python scripts/check_prose.py

sweep:
	lpcfdna sweep
	lpcfdna junction --out results/junction_report.txt

figures:
	lpcfdna figures

paper:
	$(MAKE) -C paper

clean:
	rm -rf build dist *.egg-info .pytest_cache .ruff_cache
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
	$(MAKE) -C paper clean
