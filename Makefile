all:pattern

pattern:
	rootcint -f patternDict.cpp -c pattern.h patternLinkDef.h
	g++ -c pattern.cpp -lgsl -lgslcblas `root-config --cflags --libs`
	g++ -c patternDict.cpp -lgsl -lgslcblas `root-config --cflags --libs`
	g++ pattern.o patternDict.o -o pattern.out -lgsl -lgslcblas `root-config --cflags --libs`
	g++ drawpattern.cc -o drawpattern.out `root-config --cflags --libs` 
fft:
	rootcint -f ftttestDict.cpp -c ftttest.h ftttestLinkDef.h
	g++ -I .. -fopenmp ftttest.cpp ftttestDict.cpp fftw++/fftw++.cc -o ftttest.out -lgsl -lgslcblas -lfftw3 -lfftw3_threads `root-config --cflags --libs`	
#	g++ -c ftttest.cpp -lgsl -lgslcblas -lfftw3 -lfftw3_threads `root-config --cflags --libs`
#	g++ -c ftttestDict.cpp -lgsl -lgslcblas -lfftw3 -lfftw3_threads `root-config --cflags --libs`
#	g++ ftttest.o ftttestDict.o -o ftttest.out -lgsl -lfftw3 -lfftw3_threads -lgslcblas `root-config --cflags --libs`
clean:
	rm ./pattern.out 
	rm ./drawpattern.out

