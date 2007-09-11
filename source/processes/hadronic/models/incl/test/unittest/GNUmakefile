PACKAGE = TestingFramework

# OO INCL/ABLA lib
LIBOOINCLABLA = $(G4LIB)/$(G4SYSTEM)/libG4hadronic_incl_cascade.a

#Geant4
G4HEADERS = $(G4INSTALL)/source/global/management/include

# G4 libraries:
G4LIBS = $(G4INSTALL)/lib/$(G4SYSTEM)/libG4globman.a

# CLHEP header path:
CLHEPHEADERS = /usr/local/contrib/cernsw/clhep/2.0.3.1/include/

# ROOT
ROOTINC := $(shell root-config --incdir)
ROOTLIBS := $(shell root-config --libs)

# Compiler flags
CXXFLAGS = -I. -Iinclheaders/ -I$(ROOTINC) -I./include -Wall -I../../include -Ireference/include
CXXFLAGS += -I$(G4HEADERS) -I$(CLHEPHEADERS) -DG4IOS_EXPORT

# Debug
#CXXFLAGS += -g

# Needed for CINT (dictionary generation)
CINTINCLUDES = -I. -I$(ROOTINC) -I./include

# Needed libraries
LIBS += $(ROOTLIBS)

# Headers
HDRFILES := $(shell ls include/*.hh | xargs)

# Framework infrastructure:
# HDRFILES = include/InclAblaFunctions.hh include/InclAblaTestConsoleLogWriter.hh
# HDRFILES += include/InclAblaTestDifferencePlotter.hh include/InclAblaTest.hh
# HDRFILES += include/InclAblaTestHtmlLogWriter.hh include/InclAblaTestManager.hh
# HDRFILES += include/InclAblaTestSet.hh include/InclAblaTestLogWriter.hh

# # Tests
# HDRFILES += include/TestDeutv.hh include/TestRibm.hh
# HDRFILES += include/TestDerivGaus.hh include/TestDerivMho.hh
# HDRFILES += include/TestRgaus.hh include/TestDerivWsax.hh
# HDRFILES += include/TestWsax.hh include/TestDmho.hh

LINKFILES = include/$(PACKAGE)LinkDef.h

# Obj. files
SRCOBJFILES := $(shell ls src/*.cc | xargs | sed s/\.cc/\.o/g)

# Framework infrastructure
# SRCOBJFILES = src/InclAblaTestLogWriter.o src/InclAblaTestConsoleLogWriter.o
# SRCOBJFILES += src/InclAblaTestDifferencePlotter.o src/InclAblaTest.o
# SRCOBJFILES += src/InclAblaTestHtmlLogWriter.o src/InclAblaTestManager.o
# SRCOBJFILES += src/InclAblaTestSet.o

# # Tests
# SRCOBJFILES += src/TestDeutv.o src/TestRibm.o
# SRCOBJFILES += src/TestDerivGaus.o src/TestDerivMho.o
# SRCOBJFILES += src/TestRgaus.o src/TestDerivWsax.o
# SRCOBJFILES += src/TestWsax.o src/TestDmho.o 

# ROOT dictionary
DICTOBJFILES = dict/$(PACKAGE)_dict.o

SO = lib$(PACKAGE).so
STATICLIB = lib$(PACKAGE).a

all: $(SO) RootMap 

$(SO) : $(SRCOBJFILES) dictionary
	$(CXX) -o $@ -shared $(SRCOBJFILES) $(DICTOBJFILES)  $(G4LIBS) $(LIBOOINCLABLA) $(LIBS)

$(STATICLIB) : $(SRCOBJFILES) dictionary
	$(CXX) -o $@ -static $(SRCOBJFILES) $(DICTOBJFILES) $(LIBS)

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

dictionary : $(HDRFILES) $(LINKFILES)
	mkdir -p dict
	rootcint -f dict/$(PACKAGE)_dict.cc -c $(HDRFILES) $(LINKFILES)
	$(CXX) $(CXXFLAGS) -c -o dict/$(PACKAGE)_dict.o dict/$(PACKAGE)_dict.cc

RootMap:
	rlibmap -f -r inclabla.rootmap -l $(SO) -d libCore.so -c $(LINKFILES)

# For Emacs
tags:
	etags -l c++ --declarations --members include/*.hh src/*.cc

classtree:
	ebrowse include/*.hh src/*.cc

.PHONY: clean
clean:
	rm -f $(SO)
	rm -f src/*.o
	rm -f dict/*.h
	rm -f dict/*.cc
	rm -f dict/*.o

.PHONY: cleandocs
cleandocs:
	rm -rf htmldoc/*


