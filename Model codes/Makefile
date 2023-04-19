CC = g++
CFLAGS = -O2 -Wall -std=gnu++14
LDFLAGS =
INCLUDES =

SOURCES = main.cpp
HEADERS =
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = avian.exe

avian: $(SOURCES) $(EXECUTABLE) $(HEADERS)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(INCLUDES) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean

clean:
	rm -rf *.o *.exe

