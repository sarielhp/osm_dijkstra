all: dj

LIBS= -lcairo -lpthread -lz -lexpat -lbz2

# --std=c++11 
read_file.o: read_file.cpp
	g++ -Wall -c read_file.cpp

dj: read_file.o
	g++ -o dj read_file.o $(LIBS)
