EXEC := main.out
OBJDIR := obj
SRCDIR := src
DEPDIR := $(OBJDIR)/.deps
DATADIR := data

CXX := g++
CXXFLAGS := -std=c++17 -O2 -I$(SRCDIR)
LINKFLAGS := $(CXXFLAGS)

SRCFILES := $(shell ls -A $(SRCDIR) | grep -E \.cpp$)
DEPFILES := $(SRCFILES:%.cpp=$(DEPDIR)/%.d)
OBJFILES := $(SRCFILES:%.cpp=$(OBJDIR)/%.o)

$(DEPDIR): ; @mkdir -p $@
$(DATADIR): ; @mkdir -p $@

run: $(EXEC) | $(DATADIR)
	cd $(DATADIR) && ../$(EXEC) 1.0 1e-7
	cd $(DATADIR) && ../$(EXEC) 1.0 1e-9
	cd $(DATADIR) && ../$(EXEC) 1.0 1e-11 print


clean:
	rm -rf $(OBJDIR) $(EXEC)

$(EXEC): $(OBJFILES)
	$(CXX) $(OBJDIR)/*.o -o $(EXEC) $(LINKFLAGS)

$(DEPDIR)/%.d: $(SRCDIR)/%.cpp | $(DEPDIR)
	$(CXX) -E $(CXXFLAGS) $< -MM -MT $(OBJDIR)/$(*:=.o) -MF $@
	@echo '	$(CXX) $(CXXFLAGS) -c $$(filter %.cpp,$$<) -o $$@'>>$@

$(OBJDIR)/%.o: $(DEPDIR)/%.d

include $(DEPFILES)