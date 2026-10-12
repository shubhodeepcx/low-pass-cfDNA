.PHONY: all pdf docx clean

all: pdf docx

pdf: main.pdf

main.pdf: main.tex refs.bib $(wildcard ../figures/*.png)
	pdflatex -interaction=nonstopmode -halt-on-error main.tex
	bibtex main
	pdflatex -interaction=nonstopmode -halt-on-error main.tex
	pdflatex -interaction=nonstopmode -halt-on-error main.tex

docx: main.docx

main.docx: main.tex refs.bib ieee.csl
	pandoc main.tex --bibliography=refs.bib --csl=ieee.csl --citeproc \
	  --resource-path=.:../figures -o main.docx

clean:
	rm -f *.aux *.log *.out *.bbl *.blg *.synctex.gz *.fls *.fdb_latexmk
