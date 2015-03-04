CC=g++
CFLAGS=-O3 -Wall -ffast-math -std=c++11 -fopenmp
LDFLAGS=-lm -ffast-math -std=c++11 -fopenmp
SOURCES= \
Main.cpp \
WVS.cpp

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf *.o
	rm -rf $(EXECUTABLE)
