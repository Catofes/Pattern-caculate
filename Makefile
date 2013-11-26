all:pattern

pattern:
	rootcint -f patternDict.cpp -c pattern.h patternLinkDef.h
	g++ -c pattern.cpp -lgsl -lgslcblas `root-config --cflags --libs`
	g++ -c patternDict.cpp -lgsl -lgslcblas `root-config --cflags --libs`
	g++ pattern.o patternDict.o -o pattern.out -lgsl -lgslcblas `root-config --cflags --libs`
clean:
	rm ./pattern.out 

