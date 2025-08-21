LATEX_CC=xelatex
BIBER_CC=biber

OUTPUT_DIR=build

default-target: docs

docs: main.pdf

main.pdf: main.tex preamble.tex macros.tex
	$(LATEX_CC) main $(LCC_FLAGS)
	
# And this is why we hate the LaTeX ecosystem... I'ts so.. fun.
bibliography:
	$(LATEX_CC) main
	$(BIBER_CC) main
	$(LATEX_CC) main
