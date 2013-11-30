all:pattern

pattern:
	rootcint -f patternDict.cpp -c pattern.h patternLinkDef.h
	g++ -c pattern.cpp -lgsl -lgslcblas `root-config --cflags --libs`
	g++ -c patternDict.cpp -lgsl -lgslcblas `root-config --cflags --libs`
	g++ pattern.o patternDict.o -o pattern.out -lgsl -lgslcblas `root-config --cflags --libs`
	g++ drawpattern.cc -o drawpattern.out `root-config --cflags --libs` 
	rootcint -f ftttestDict.cpp -c ftttest.h ftttestLinkDef.h
	g++ -c ftttest.cpp -lgsl -lgslcblas `root-config --cflags --libs`
	g++ -c ftttestDict.cpp -lgsl -lgslcblas `root-config --cflags --libs`
	g++ ftttest.o ftttestDict.o -o ftttest.out -lgsl -lgslcblas `root-config --cflags --libs`
	rm pattern.o patternDict.o patternDict.cpp patternDict.h ftttest.o ftttestDict.o ftttestDict.cpp ftttestDict.h
clean:
	rm ./pattern.out 
	rm ./drawpattern.out

