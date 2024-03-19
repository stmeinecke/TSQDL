all:
	g++ -std=c++11 TSQDL.cpp -o TSQDL -Ofast -I $$PWD/incl/ -lfftw3 -lm
