src = $(wildcard src/*.cpp)
obj = $(src:.cpp=.o)
tex = $(proposal:.tex=.aux)
PAPER=$(wildcard proposal/*.tex)  # set the path to your TeX file here
CC = g++

CXXFLAGS = -g -std=c++11

main: $(obj)
	$(CC) -o $@ $^ $(CXXFLAGS)

.PHONY: clean

clean:
	rm -f $(obj) main
	rm -f $(tex) prop

prop: $(tex)
	pdflatex $(PAPER)
	bibtex $(tex)