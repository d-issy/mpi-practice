all: clear matrix

.PHONY: matrix clear
matrix:
	mpic++ matrix.cpp -o $@

clear:
	rm -rf matrix
