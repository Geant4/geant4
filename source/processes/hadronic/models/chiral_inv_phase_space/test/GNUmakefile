# ----------------------------------------------------------------
# File: CHIPStest makefile
# Author: Mikhail Kossov (based on Geant4 examples & CMSExample of I.G.)
# Date: 26/08/99
# Modifications:
# ----------------------------------------------------------------

#
# Target selection
#
name := CHIPStest
ifndef G4TARGET
  G4TARGET := $(name)
endif

OSPACE_DIR = /afs/cern.ch/sw/lhcxx/specific/sun/ObjectSpace/pro
STL_DIR = ${OSPACE_DIR}/ospace/stl
STD_DIR = ${OSPACE_DIR}/ospace/std
CC_DEFINES = -DOS_SOLARIS_2_5 -DOS_NEW_CHECK -DOS_STL_ASSERT -DOS_NO_WSTRING  -DOS_NO_ALLOCATORS
CC_DEFINES+= -DHEP_HAVE_STL 
#CC_DEFINES+= -DCLHEP_MAX_MIN_DEFINED 
CC_DEFINES+= -DOS_STL_FUNCTION_H -DOS_OMIT_MIN_MAX

#
# Standard G4 things
#
MAKEFLAGS= --no-print-directory

include $(G4INSTALL)/config/binmake.gmk
CPPFLAGS +=  -I$(OSPACE_DIR) -I$(STL_DIR) -I$(STD_DIR) $(CC_DEFINES)

# So that HP-UX/CC find My*, G4Allocator things. 
CXXFLAGS  += -g

EXEDIR := $(G4TMPDIR)/exe


#
# CHIPStest include structure
#
CPPFLAGS += -I. -ICUtils -ICData -ICBody -ICNucleus -ICG4Interface 


#
#Add paths & libraries for the Persistentcy
#
ifndef TEST_NO_OBJY
CPPFLAGS += -I$(OBJY_DIR)/include -I$(HEP_ODBMS_DIR)/include \
            -DHep_HAS_BUILDIN_BOOL_CONSTANTS

LDFLAGS  += -L$(OBJY_DIR)/lib -L$(HEP_ODBMS_DIR)/lib
LDLIBS   += -loo.4.02 -lHepODBMS
CPPFLAGS += -IPersistency -DTEST_USE_OBJY
endif


#
# Define source files for the gmake dependencies
#
sources1 := $(wildcard $(G4TARGET).cc)
sources2 := $(wildcard CUtils/*.cc)
sources3 := $(wildcard CData/*.cc)
sources4 := $(wildcard CBody/*.cc)
sources5 := $(wildcard CNucleus/*.cc)
sources6 := $(wildcard CG4Interface/*.cc)
ifndef TEST_NO_OBJY
  sources7 := $(wildcard Persistency/*.cc)
endif

ifeq ($(G4TARGET),CHIPStest)
  sources  := $(sources1) $(sources2) $(sources3) $(sources4) $(sources5) $(sources6) $(sources7)
endif
#ifeq ($(G4TARGET),FillDB)
#  sources := $(sources1) $(sources3) $(sources6) $(sources7)
#endif

javasources := $(wildcard GUI/*.java)

nobjyobj := $(patsubst %.cc,$(EXEDIR)/%.o,$(sources))
javaobjects := $(patsubst %.java,%.class,$(javasources))

objects  := $(nobjyobj) 

nobjydep     := $(patsubst %.o,%.d,$(objects))
dependencies := $(nobjydep)

ifndef TEST_NO_OBJY
  objysources :=$(wildcard Persistency/*.C)
  objyobj  := $(patsubst %.C,$(EXEDIR)/%.o,$(objysources))
  objects  += $(objyobj)
  objydep      := $(patsubst %.o,%.d,$(objyobj))
  dependencies += $(objydep)
endif


# Make $(G4TARGET) executable.
$(G4BINDIR)/$(G4TARGET): $(objects) $(javaobjects)

java: $(javaobjects)

%.class : %.java
	@echo Compiling $(<F)
	@cd GUI; javac $(<F)

# Make the .o files(s).
$(EXEDIR)/%.o : %.cc
ifdef CPPVERBOSE
	@echo Compiling $< ...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(APPFLAGS) \
	-c -o $@ $<
else
	@echo Compiling $< ...
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(APPFLAGS) \
	-c -o $@ $<
endif

$(EXEDIR)/%.o : %.C
ifdef CPPVERBOSE
	@echo Compiling $< ...
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(APPFLAGS) \
	-c -o $@ $<
else
	@echo Compiling $< ...
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(APPFLAGS) \
	-c -o $@ $<
endif

# Make the .d file(s) and include it(them).
$(EXEDIR)/%.d: %.cc 
	@if [ ! -d $(G4TMP)/$(G4SYSTEM) ] ; then mkdir $(G4TMP)/$(G4SYSTEM) ;fi
	@if [ ! -d $(G4TMPDIR) ] ; then mkdir $(G4TMPDIR) ;fi
	@if [ ! -d $(G4TMPDIR)/exe ] ; then mkdir $(G4TMPDIR)/exe ;fi
	@if [ ! -d $(EXEDIR)/CUtils ] ; then mkdir $(EXEDIR)/CUtils ;fi
	@if [ ! -d $(EXEDIR)/CData ] ; then mkdir $(EXEDIR)/CData ;fi
	@if [ ! -d $(EXEDIR)/CBody ] ; then mkdir $(EXEDIR)/CBody ;fi
	@if [ ! -d $(EXEDIR)/CNucleus ] ; then mkdir $(EXEDIR)/CNucleus ;fi
	@if [ ! -d $(EXEDIR)/CG4Interface ] ; then mkdir $(EXEDIR)/CG4Interface ;fi
ifndef TEST_NO_OBJY
	@if [ ! -d $(EXEDIR)/Persistency ] ; then mkdir $(EXEDIR)/Persistency ;fi
endif
	@echo Making dependency file $@...
ifeq ($(G4SYSTEM),G++-Linux)
	@($(ECHO) -n $(G4TMPDIR)/exe/ ; \
	g++ -MM $(RWINCGPP) $(CPPFLAGS) $< ) | sed 's!$(G4TMPDIR)/exe/$*.o!& $@!' >$@
else
	@($(ECHO) $(G4TMPDIR)/\\c ; \
	g++ -MM $(RWINCGPP) $(CPPFLAGS) $< ) | sed 's!$(G4TMPDIR)/exe/$*.o!& $@!' >$@
endif
$(EXEDIR)/%.d: %.C
	@if [ ! -d $(G4TMP)/$(G4SYSTEM) ] ; then mkdir $(G4TMP)/$(G4SYSTEM) ;fi
	@if [ ! -d $(G4TMPDIR) ] ; then mkdir $(G4TMPDIR) ;fi
	@if [ ! -d $(G4TMPDIR)/exe ] ; then mkdir $(G4TMPDIR)/exe ;fi
	@if [ ! -d $(EXEDIR)/CUtils ] ; then mkdir $(EXEDIR)/CUtils ;fi
	@if [ ! -d $(EXEDIR)/CData ] ; then mkdir $(EXEDIR)/CData ;fi
	@if [ ! -d $(EXEDIR)/CBody ] ; then mkdir $(EXEDIR)/CBody ;fi
	@if [ ! -d $(EXEDIR)/CNucleus ] ; then mkdir $(EXEDIR)/CNucleus ;fi
	@if [ ! -d $(EXEDIR)/CG4Interface ] ; then mkdir $(EXEDIR)/CG4Interface ;fi
ifndef TEST_NO_OBJY
	@if [ ! -d $(EXEDIR)/Persistency ] ; then mkdir $(EXEDIR)/Persistency ;fi
endif
	@echo Making dependency file $@...
ifeq ($(G4SYSTEM),G++-Linux)
	@($(ECHO) -n $(G4TMPDIR)/exe/ ; \
	g++ -MM $(RWINCGPP) $(CPPFLAGS) $< ) | sed 's!$(G4TMPDIR)/exe/$*.o!& $@!' >$@
else
	@($(ECHO) $(G4TMPDIR)/\\c ; \
	g++ -MM $(RWINCGPP) $(CPPFLAGS) $< ) | sed 's!$(G4TMPDIR)/exe/$*.o!& $@!' >$@
endif
-include $(dependencies)

