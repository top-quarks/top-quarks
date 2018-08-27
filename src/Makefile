cern: main.cpp Makefile polarmodule.hpp input.hpp io.hpp score.hpp old_pairs.hpp load.hpp investigate.hpp density.hpp new_pairs.hpp
	g++ main.cpp -O3 -g -std=c++11 -o cern -Wfatal-errors

run: SHELL:=/bin/bash
run: cern
	time ./cern
