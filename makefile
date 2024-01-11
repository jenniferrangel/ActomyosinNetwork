# **************************************************
# Variables to control Makefile operation
CXX= -g++ -fopenmp -static-libstdc++
CXXFLAGS= -Wall -g -O3
#*****************************************************
# Targets needed to bring the executable up to date
all: program 

program:folder main.o coord.o node.o filament.o network.o
	$(CXX) $(CXXFLAGS) -o program main.o coord.o node.o filament.o network.o
folder: 
		mkdir -p ./DataOutput ./Animation
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

coord.o: coord.cpp
	$(CXX) $(CXXFLAGS) -c coord.cpp

node.o: node.cpp
	$(CXX) $(CXXFLAGS) -c node.cpp

filament.o: filament.cpp
	$(CXX) $(CXXFLAGS) -c filament.cpp

network.o: network.cpp
	$(CXX) $(CXXFLAGS) -c network.cpp

clean: wipe
		rm -rf *o program

wipe: