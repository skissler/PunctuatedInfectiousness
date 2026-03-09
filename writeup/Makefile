all: main.pdf

main.pdf: main.tex introduction.tex model.tex uncontrolled.tex interventions.tex identifiability.tex discussion.tex acknowledgements.tex supplementary.tex references.bib
	pdflatex -interaction=nonstopmode main.tex
	bibtex main
	pdflatex -interaction=nonstopmode main.tex
	pdflatex -interaction=nonstopmode main.tex

clean:
	rm -f main.aux main.bbl main.blg main.log main.out main.pdf main.run.xml

.PHONY: all clean
