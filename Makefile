all:pattern

pattern:
	g++ -o pattern.out pattern.cpp -lgsl -lgslcblas 
clean:
	rm ./pattern.out 

