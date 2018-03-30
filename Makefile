CC := g++
CFLAGS := -std=c++17 -O3 -fopenmp -pedantic -Wall -m64

PROJECTNAME := AstroQUT

SOURCEDIR := src
SOURCES := $(shell find $(SOURCEDIR) -name '*.cpp')
HEADERDIR := include
EXTERNALHEADERDIR := extern_include/eigen-eigen-5a0156e40feb
OBJECTDIR := obj
OBJECTS := $(addprefix $(OBJECTDIR)/,$(SOURCES:%.cpp=%.o))
DEPENDS := $(addprefix $(OBJECTDIR)/,$(SOURCES:%.cpp=%.d))


all: $(PROJECTNAME)

rebuild: clean $(PROJECTNAME)

$(PROJECTNAME): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -fopenmp -o $(PROJECTNAME)

-include $(DEPENDS)

$(OBJECTDIR)/%.o: %.cpp
	mkdir -p $(OBJECTDIR)/$(dir $<) && $(CC) $(CFLAGS) -I $(HEADERDIR) -I $(EXTERNALHEADERDIR) -c $< -o $@

.PHONY: clean

clean:
	rm -rf $(OBJECTDIR) $(PROJECTNAME)
