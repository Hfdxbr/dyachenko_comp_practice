EXEC := main.out
OBJDIR := obj
SRCDIR := src
DEPDIR := $(OBJDIR)/.deps
DATADIR := data

CXX := g++
CXXFLAGS := -std=c++17 -I$(SRCDIR)
LINKFLAGS := $(CXXFLAGS)

SRCFILES := $(shell ls -A $(SRCDIR) | grep -E \.cpp$)
DEPFILES := $(SRCFILES:%.cpp=$(DEPDIR)/%.d)
OBJFILES := $(SRCFILES:%.cpp=$(OBJDIR)/%.o)

$(DEPDIR): ; @mkdir -p $@
$(DATADIR): ; @mkdir -p $@

run: $(EXEC) | $(DATADIR)
	./$(EXEC)

clean:
	rm -rf $(OBJDIR) $(EXEC)

$(EXEC): $(OBJFILES)
	$(CXX) $(OBJDIR)/*.o -o $(EXEC) $(LINKFLAGS)

$(DEPDIR)/%.d: $(SRCDIR)/%.cpp | $(DEPDIR)
	$(CXX) -E $(CXXFLAGS) $< -MM -MT $(OBJDIR)/$(*:=.o) -MF $@
	@echo '	$(CXX) $(CXXFLAGS) -c $$(filter %.cpp,$$<) -o $$@'>>$@

$(OBJDIR)/%.o: $(DEPDIR)/%.d

include $(DEPFILES)