# $Id: GNUmakefile,v 1.1 1999-01-08 16:33:56 gunter Exp $

MAKEFLAGS= --no-print-directory

ifndef TESTTARGET
 TESTTARGET := test19
endif

ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif
 
ifndef G4INSTALL
  G4INSTALL = ../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

G4TMPDIR := $(G4TMP)/$(G4SYSTEM)/$(G4TARGET)
G4BINDIR := $(G4BIN)/$(G4SYSTEM)

ifdef G4VIS_DEBUG
  G4VIS_OPTIMISESTRING := G4DEBUG=$(G4VIS_DEBUG)
endif

include $(G4INSTALL)/config/G4UI_USE.gmk
include $(G4INSTALL)/config/G4VIS_USE.gmk
include $(G4INSTALL)/config/interactivity.gmk

# Override some variables for binmake.gmk.
#
ifeq ($(G4TARGET), BaBarCallListTest)
  CPPFLAGS += -Ibabar_utils/include
  CPPFLAGS += -I.
endif
CPPFLAGS += -Itest_utils/include

CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/run/include
CPPFLAGS += -I$(G4BASE)/event/include
CPPFLAGS += -I$(G4BASE)/tracking/include
CPPFLAGS += -I$(G4BASE)/track/include
CPPFLAGS += -I$(G4BASE)/global/HEPGeometry/include
CPPFLAGS += -I$(G4BASE)/geometry/volumes/include
CPPFLAGS += -I$(G4BASE)/geometry/management/include
CPPFLAGS += -I$(G4BASE)/particles/management/include
CPPFLAGS += -I$(G4BASE)/processes/management/include
CPPFLAGS += -I$(G4BASE)/processes/electromagnetic/utils/include
CPPFLAGS += -I$(G4BASE)/global/HEPRandom/include
CPPFLAGS += -I$(G4BASE)/materials/include
CPPFLAGS += -I$(G4BASE)/digits+hits/detector/include
CPPFLAGS += -I$(G4BASE)/digits+hits/hits/include
CPPFLAGS += -I$(G4BASE)/digits+hits/digits/include
CPPFLAGS += -I$(G4BASE)/intercoms/include
CPPFLAGS += -I$(G4BASE)/graphics_reps/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/CSG/include

ifdef G4VIS_USE_OPACS
  CPPFLAGS += -I$(COROOT)/include -I$(HOROOT)/include
endif

LDFLAGS   += -L$(G4LIB)/$(G4SYSTEM)
LDFLAGS   += -L$(CLHEP_LIB_DIR)
LDFLAGS   += -L$(RWLIBDIR)

# test19/test_utils...
ifeq ($(G4TARGET), BaBarCallListTest)
  LOADLIBS := -lG4babar_utils
endif
LOADLIBS += -lG4vis_utils

ifdef G4ATLAS
  LOADLIBS += -lG4ATLAS
endif

LOADLIBS += $(VISLIBS) $(UILIBS)
LOADLIBS += -lG4run -lG4event
LOADLIBS += -lG4tracking
LOADLIBS += -lG4particles
LOADLIBS += -lG4processes
LOADLIBS += -lG4track
LOADLIBS += -lG4materials -lG4digits+hits -lG4geometry
LOADLIBS += -lG4graphics_reps -lG4intercoms  -lG4global
LOADLIBS += $(INTYLIBS)
LOADLIBS += -l$(CLHEP_LIB) -l$(RWLIB) -lm

sources := $(wildcard $(G4TARGET).cc)
objects := $(patsubst %.cc,$(G4TMPDIR)/%.o,$(sources))
dependencies := $(patsubst %.cc,$(G4TMPDIR)/%.d,$(sources))

libs := $(patsubst -l%,lib%.a,$(filter -l%, $(LOADLIBS) $(LDFLAGS)))
libpath := $(patsubst -L%,%,$(filter -L%,$(LOADLIBS) $(LDFLAGS)))

## --------------------- vpaths ----------------------------------------------

vpath %.hh  $(patsubst -I%,%,$(filter -I%,$(CPPFLAGS)))
vpath %.a   $(libpath)

### --------------------- vpaths ----------------------------------------------

# Make $(G4TARGET) executable.
#$(G4SYSTEM)/$(G4TARGET): $(objects) $(libs)
$(G4SYSTEM)/$(G4TARGET): $(objects)
ifndef link
ifndef vis
	@cd $(G4BASE)/global &&          $(MAKE) global
	@cd $(G4BASE)/intercoms &&       $(MAKE)
	@cd $(G4BASE)/interfaces &&      $(MAKE)
	@cd $(G4BASE)/event &&           $(MAKE)
	@cd $(G4BASE)/digits+hits &&     $(MAKE) global
	@cd $(G4BASE)/tracking &&        $(MAKE)
	@cd $(G4BASE)/track &&           $(MAKE)
	@cd $(G4BASE)/particles &&       $(MAKE) global
	@cd $(G4BASE)/materials &&       $(MAKE) 
	@cd $(G4BASE)/processes &&       $(MAKE) global
	@cd $(G4BASE)/run &&             $(MAKE)
	@cd $(G4BASE)/geometry &&        $(MAKE) global
ifdef G4ATLAS
	@cd $(G4ATLAS) &&        $(MAKE)
endif
endif
	@cd $(G4BASE)/graphics_reps &&   $(MAKE) $(G4VIS_OPTIMISESTRING)
	@cd $(G4BASE)/visualization &&   $(MAKE) $(G4VIS_OPTIMISESTRING)
ifeq ($(G4TARGET), BaBarCallListTest)
	@cd babar_utils &&           $(MAKE) $(G4VIS_OPTIMISESTRING)
endif
	@cd test_utils &&            $(MAKE) $(G4VIS_OPTIMISESTRING)
endif
	@if [ ! -d $(G4BINDIR) ] ; then mkdir $(G4BINDIR) ;fi
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) \
	 -o $(G4BINDIR)/$(G4TARGET) $(objects) $(LDFLAGS) $(LOADLIBS)

# Make the .o files(s).
$(G4TMPDIR)/%.o : %.cc
	@rm -f $(G4TMPDIR)/%.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(APPFLAGS) \
	 -c -o $(G4TMPDIR)/$(G4TARGET).o $(G4TARGET).cc

# Make the .d file(s) and include it(them).
$(G4TMPDIR)/%.d: %.cc
	@if [ ! -d $(G4TMP)/$(G4SYSTEM) ] ; then mkdir $(G4TMP)/$(G4SYSTEM) ;fi
	@if [ ! -d $(G4TMPDIR) ] ; then mkdir $(G4TMPDIR) ;fi
	@echo Making dependency for file $<...
ifeq ($(G4SYSTEM),Linux-g++)
	@echo -n $(G4TMPDIR)/ >$(G4TMP)/$(USER)_$(G4SYSTEM).d
else
	@echo $(G4TMPDIR)/\\c >$(G4TMP)/$(USER)_$(G4SYSTEM).d
endif
	@g++ -MM $(CPPFLAGS) $< >>$(G4TMP)/$(USER)_$(G4SYSTEM).d
	@sed 's!$(G4TMPDIR)/$*.o!& $@!' $(G4TMP)/$(USER)_$(G4SYSTEM).d >$@
	@rm $(G4TMP)/$(USER)_$(G4SYSTEM).d
-include $(dependencies)

.PHONY: quick anew clean clean_libs clean_all clean_all_tildas

quick:
	$(MAKE) anew vis=1

anew:
	-rm -f $(G4BINDIR)/$(G4TARGET)
	$(MAKE)

clean:
	@echo Cleaning up ...
	-rm -f FR*
	-rm -f g4*
	-rm -f *.list
	@rm -f $(G4BINDIR)/$(G4TARGET)
ifdef CFRONT
	@rm -rf $(G4TREP)/exec
endif

clean_libs:
	@cd $(G4BASE)/global; $(MAKE) clean_libs
	@cd $(G4BASE)/intercoms; $(MAKE) clean_libs
	@cd $(G4BASE)/interfaces; $(MAKE) clean_libs
	@cd $(G4BASE)/particle+matter; $(MAKE) clean_libs
	@cd $(G4BASE)/tracking; $(MAKE) clean_libs
	@cd $(G4BASE)/track; $(MAKE) clean_libs
	@cd $(G4BASE)/event; $(MAKE) clean_libs
	@cd $(G4BASE)/geometry; $(MAKE) clean_libs
	@cd $(G4BASE)/digits+hits; $(MAKE) clean_libs
	@cd $(G4BASE)/graphics_reps; $(MAKE) clean_libs
ifdef G4ATLAS
	@cd $(G4ATLAS); $(MAKE) clean_libs
endif
	@cd ..;  $(MAKE) clean_libs
	@cd utils;  $(MAKE) clean_libs
	@$(MAKE) clean

clean_all:
	@cd $(G4BASE)/global; $(MAKE) clean
	@cd $(G4BASE)/intercoms; $(MAKE) clean
	@cd $(G4BASE)/interfaces; $(MAKE) clean
	@cd $(G4BASE)/particle+matter; $(MAKE) clean
	@cd $(G4BASE)/tracking; $(MAKE) clean
	@cd $(G4BASE)/track; $(MAKE) clean
	@cd $(G4BASE)/event; $(MAKE) clean
	@cd $(G4BASE)/geometry; $(MAKE) clean
	@cd $(G4BASE)/digits+hits; $(MAKE) clean
	@cd $(G4BASE)/graphics_reps; $(MAKE) clean
ifdef G4ATLAS
	@cd $(G4ATLAS); $(MAKE) clean
endif
	@cd ..; $(MAKE) clean
	@cd utils; $(MAKE) clean
	@$(MAKE) clean
ifdef CFRONT
	@rm -rf $(G4TREP)
endif

clean_all_tildas:
	find $(G4BASE)/../../ -name '*~' -exec rm -i {} \;
