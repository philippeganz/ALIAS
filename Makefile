CC := mpic++
CFLAGS := -std=c++1z -O2 -fopenmp -pedantic -Wall
DFLAGS := -MM -std=c++1z

PROJECTNAME := ASTROQUT

SOURCEDIR := src
SOURCES := $(shell find $(SOURCEDIR) -name '*.cpp')
HEADERDIR := include
OBJECTDIR := obj
OBJECTS := $(addprefix $(OBJECTDIR)/,$(SOURCES:%.cpp=%.o))
DEPENDS := $(addprefix $(OBJECTDIR)/,$(SOURCES:%.cpp=%.d))


all: $(PROJECTNAME)

rebuild: clean $(PROJECTNAME)

$(PROJECTNAME): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -fopenmp -o $(PROJECTNAME)

-include $(DEPENDS)

$(OBJECTDIR)/%.o: %.cpp
	mkdir -p $(OBJECTDIR)/$(dir $<) && $(CC) $(CFLAGS) -I $(HEADERDIR) -c $< -o $@
	$(CC) $(DFLAGS) -I $(HEADERDIR) $< > $(OBJECTDIR)/$(*D)/$(*F).d
	@mv -f $(OBJECTDIR)/$(*D)/$(*F).d $(OBJECTDIR)/$(*D)/$(*F).d.tmp
	@sed -e 's|.*:|$(OBJECTDIR)/$(*D)/$(*F).o:|' < $(OBJECTDIR)/$(*D)/$(*F).d.tmp > $(OBJECTDIR)/$(*D)/$(*F).d
	@sed -e 's/.*://' -e 's/\\$$//' < $(OBJECTDIR)/$(*D)/$(*F).d.tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(OBJECTDIR)/$(*D)/$(*F).d
	@rm -f $(OBJECTDIR)/$(*D)/$(*F).d.tmp

.PHONY: clean

clean:
rm -rf $(OBJECTDIR) $(PROJECTNAME)
