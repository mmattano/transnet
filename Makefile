# TransNet: build the networks, run the notebooks, build the docs.
#
# The networks need internet and BRENDA credentials (BRENDA_EMAIL,
# BRENDA_PASSWORD) and take 30-90 minutes each; they checkpoint, so an
# interrupted build resumes. Everything else is fast and offline.

PYTHON ?= python3
ORGANISMS ?= mouse human rat yeast ecoli
# The teaching order, which the filenames deliberately no longer carry. The
# same order is in docs/source/index.rst, and tests/test_notebooks.py fails if
# the two disagree. external_annotation needs the network, so `make notebooks`
# leaves it out; run it by hand.
WALKTHROUGHS := build_network responsive_network reaction_regulation \
                regulatory_paths temporal_and_hubs compare_conditions \
                network_topology export_network transcription_factors
STUDIES := $(wildcard notebooks/studies/*.py)

.PHONY: help networks notebooks studies docs test test-all clean

help:
	@echo "make networks   rebuild data/<organism>/latest by hand for: $(ORGANISMS)"
	@echo "make notebooks  run the offline walkthroughs"
	@echo "make studies    run the studies (needs the networks)"
	@echo "make docs       build the documentation into docs/build/html"
	@echo "make test       run the test suite (offline, fast)"
	@echo "make test-all   also run the notebooks end to end"

networks:
	@for organism in $(ORGANISMS); do \
		echo "== $$organism"; \
		$(PYTHON) maintenance/build_networks.py --organisms $$organism --brenda || exit 1; \
	done

notebooks:
	@for name in $(WALKTHROUGHS); do \
		echo "== $$name"; \
		MPLBACKEND=Agg $(PYTHON) notebooks/walkthroughs/$$name.py > /dev/null || exit 1; \
	done
	@echo "all walkthroughs ran"

# The studies are executed into docs/source/studies/*.ipynb with their outputs,
# which is what the documentation renders -- they need built networks and take
# minutes, so the docs build does not re-run them.
studies:
	@mkdir -p docs/source/studies
	@for study in $(STUDIES); do \
		name=$$(basename $$study .py); \
		echo "== $$name"; \
		(cd notebooks/studies && jupytext --to ipynb --set-kernel transnet \
			--execute $$name.py -o ../../docs/source/studies/$$name.ipynb) || exit 1; \
	done
	$(PYTHON) maintenance/refresh_doc_figures.py

docs:
	$(MAKE) -C docs html

test:
	$(PYTHON) -m pytest tests -q -m "not network and not slow"
	$(PYTHON) maintenance/verify_references.py --offline-ok

test-all:
	$(PYTHON) -m pytest tests -q -m "not network"

clean:
	rm -rf docs/build notebooks/exports notebooks/walkthroughs/exports
	find . -name __pycache__ -not -path "./venv/*" -exec rm -rf {} +
