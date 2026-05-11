.PHONY: reproduce paper clean check-no-archives arxiv-source

reproduce:
	python scripts/reproduce_all.py --out results

paper:
	cd paper && pdflatex main.tex && biber main && pdflatex main.tex && pdflatex main.tex

check-no-archives:
	@if find . \( -name '*.zip' -o -name '*.tar' -o -name '*.tar.gz' \) | grep -q .; then echo 'Archives found in repository tree:'; find . \( -name '*.zip' -o -name '*.tar' -o -name '*.tar.gz' \); exit 1; else echo 'OK: no archives inside repository tree'; fi

arxiv-source:
	rm -rf build/arxiv_source build/arxiv_source.zip
	mkdir -p build/arxiv_source/figures
	cp paper/main.tex build/arxiv_source/main.tex
	cp paper/refs.bib build/arxiv_source/refs.bib
	cp paper/figures/*.png build/arxiv_source/figures/
	cd build/arxiv_source && zip -r ../arxiv_source.zip .

clean:
	find . \( -name "*.aux" -o -name "*.log" -o -name "*.out" -o -name "*.toc" -o -name "*.bcf" -o -name "*.run.xml" -o -name "*.blg" -o -name "*.bbl" -o -name "*.fls" -o -name "*.fdb_latexmk" \) -delete
	rm -rf results results_minimal build
