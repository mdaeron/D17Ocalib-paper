default: all

all: cleanup linear process alpha18 DD17Ovsklna18 SLAP sensitivity labs export-toml ftir

export-toml:
	cp output/*toml manuscript/results/

cleanup:
	rm -rf output/*
	rm -rf manuscript/results/*
	rm -rf manuscript/plots/*

edist:
	-open -a Edist manuscript/D17Ocalib-manuscript.typ

typst:
	-open -a TeXShop manuscript/D17Ocalib-manuscript.pdf
	-zed manuscript/D17Ocalib-manuscript.typ
	-cd manuscript; typst watch --ignore-system-fonts --font-path ./fonts D17Ocalib-manuscript.typ

linear:
	uv run code/linear-conversion.py
	cp output/linear-conversion.pdf  manuscript/plots/
	cp output/linear-conversion.toml  manuscript/results/

harvest:
	uv run code/harvest_rawdata.py

process:
	uv run code/process_analyses.py
	cp output/DVHW.pdf  manuscript/plots/
	cp output/*.toml  manuscript/results/

alpha18:
	uv run code/alpha18.py
	cp output/alpha18.pdf  manuscript/plots/

DD17Ovsklna18:
	uv run code/DD17O_vs_klna18.py
	cp output/DD17O_vs_klna18*.pdf  manuscript/plots/
	cp output/*.toml  manuscript/results/

SLAP:
	uv run code/SLAP_adjustments.py
	cp output/*.toml  manuscript/results/

sensitivity:
	uv run code/sensitivity.py
	cp output/sensitivity.pdf  manuscript/plots/
	cp output/*.toml  manuscript/results/

labs:
	uv run code/labs.py
	cp output/labs.pdf  manuscript/plots/
	cp output/interlab-calibrations.pdf  manuscript/plots/
	cp output/*.toml  manuscript/results/

ftir:
	uv run code/ftir.py
	cp output/ftir.pdf  manuscript/plots/

export-zip:
	rm -rf ../D17Ocalib-paper.zip
	git archive --format=zip --prefix=D17Ocalib-paper/ --output=../D17Ocalib-paper.zip HEAD -- \
		code \
		inputs \
		justfile \
		manuscript/authors.toml \
		manuscript/D17Ocalib-manuscript.pdf \
		manuscript/D17Ocalib-manuscript.typ \
		manuscript/fonts manuscript/gpl.csl \
		manuscript/images/SEM-lite.pdf \
		manuscript/plots \
		manuscript/results \
		manuscript/refs.bib  \
		manuscript/template.typ \
		output \
		pyproject.toml \
		LICENSE \
		readme.md \
		uv.lock

export-dir:
	rm -rf ../D17Ocalib-paper
	mkdir ../D17Ocalib-paper
	git archive HEAD -- \
		code \
		inputs \
		justfile \
		manuscript/plots \
		manuscript/results \
		output \
		pyproject.toml \
		LICENSE \
		readme.md \
		uv.lock \
		| tar -x -C ../D17Ocalib-paper
