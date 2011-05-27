# $Id: GNUmakefile,v 1.1 1999-01-08 16:33:38 gunter Exp $
# ----------------------------------------------------------------
# Makes test program in environment variable TESTTARGET.
# ----------------------------------------------------------------
# TESTTARGET is created according to G4SYSTEM. G.Cosmo, 11/6/96.

TESTTARGET := FitMC
include ../../../../../../global/GNUmakefile.architecture

# Override some variables from GNUmakefile.architecture.
INCFLAGS := -I../../include \
            -I$(CLHEP_INCLUDE_DIR) \
            -I$(RWINC)
             
CPPFLAGS := $(INCFLAGS)

ifdef G4_HAVE_BOOL
  CPPFLAGS += -DG4_HAVE_BOOL
endif
 
LDFLAGS  += -L../../$(G4SYSTEM) \
            -L$(CLHEP_LIB_DIR) \
            -L$(RWLIBDIR)

LOADLIBS := \
            -lG4hadronic_HE \
            -lCLHEP \
            -lrwtool  -lm 
            

possibletargets := $(basename $(wildcard *.cc))

sources := $(wildcard $(TESTTARGET).cc)
objects := $(patsubst %.cc,$(G4SYSTEM)/%.o,$(sources))
dependencies := $(patsubst %.cc,$(G4SYSTEM)/%.d,$(sources))

ifndef G4TESTEXECDIR
  G4TESTEXECDIR := $(G4SYSTEM)
endif
  
# Make $(TESTTARGET) executable.
$(G4TESTEXECDIR)/$(TESTTARGET): $(objects)
	@cd ../../../../../../processes/hadronic/models/high_energy;  $(MAKE)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(G4TESTEXECDIR)/$(TESTTARGET) $(objects) $(LDFLAGS) $(LOADLIBS)
	@echo >$(G4SYSTEM)/$(TESTTARGET).cxxflags $(CXXFLAGS)

# Make the .d file(s) and include it(them).
$(G4SYSTEM)/%.d: %.cc
	@if [ ! -d $(G4SYSTEM) ] ; then mkdir $(G4SYSTEM) ;fi
	@echo $(LOCATION)
	@echo Making dependency file $@...
	@echo $(G4SYSTEM)/\\c >/tmp/$(USER)_$(G4SYSTEM).d
	@g++ -MM $(CPPFLAGS) $< >/tmp/$(USER)_$(G4SYSTEM).d
	@sed 's!$(G4SYSTEM)/$*.o!& $@!' /tmp/$(USER)_$(G4SYSTEM).d >$@
	@rm /tmp/$(USER)_$(G4SYSTEM).d
-include $(dependencies)

# Make the .o files(s).
$(G4SYSTEM)/%.o : %.cc
	@rm -f $(G4SYSTEM)/*.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(APPFLAGS) -c -o $(G4SYSTEM)/$(*F).o $(TESTTARGET).cc

clean:
	-rm -f $(G4SYSTEM)/$(possibletargets)
	-rm -f $(G4SYSTEM)/core
	-rm -f $(G4SYSTEM)/*~ $(G4SYSTEM)/*.o $(G4SYSTEM)/*.d
	-rm -f $(G4SYSTEM)/*.cxxflags
#	@rm -r -f $(wastedirs)
#	@rm -rf $(G4SYSTEM)







