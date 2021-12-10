src = $(wildcard src/*.cc)
obj = $(src:.cc=.o)

CC = mpic++

FLAGS = -g -W -Wall -std=c++14 -fopenmp -DPRINT

PROPOSAL = proposal

main: $(obj)
	$(CC) $(FLAGS) -o main $(obj)

%.o : %.cc
	$(CC) -c $(FLAGS) $< -o $@

.PHONY: clean

clean:
	rm -f $(obj) main
	rm -f *.aux *.lof *.log *.lot *.toc *.bbl *.blg *.pdf

run: $(PROPOSAL).pdf

$(PROPOSAL).pdf: $(PROPOSAL).bbl $(PROPOSAL).tex
	pdflatex $(PROPOSAL).tex -draftmode
	pdflatex $(PROPOSAL).tex

$(PROPOSAL).bbl: $(PROPOSAL).aux
	bibtex $(PROPOSAL).aux

$(PROPOSAL).aux: $(PROPOSAL).bib
	pdflatex $(PROPOSAL).tex -draftmode
	pdflatex $(PROPOSAL).tex -draftmode
