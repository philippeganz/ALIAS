CC := g++
CFLAGS := -std=c++17 -O3 -fopenmp -pedantic -Wall -m64 -march=native -mno-sse5

PROJECTNAME := ALIAS

SOURCEDIR := src
SOURCES := $(shell find $(SOURCEDIR) -name '*.cpp')
HEADERDIR := include
OBJECTDIR := obj
OBJECTS := $(addprefix $(OBJECTDIR)/,$(SOURCES:%.cpp=%.o))

all: $(PROJECTNAME)

rebuild: clean $(PROJECTNAME)

$(PROJECTNAME): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -fopenmp -o $(PROJECTNAME)

$(OBJECTDIR)/%.o: %.cpp
	mkdir -p $(OBJECTDIR)/$(dir $<) && $(CC) $(CFLAGS) -I $(HEADERDIR) -c $< -o $@

.PHONY: clean

clean:
	rm -rf $(OBJECTDIR) $(PROJECTNAME)
