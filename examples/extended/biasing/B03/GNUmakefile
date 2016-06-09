# $Id: GNUmakefile,v 1.11 2003/06/30 16:40:27 dressel Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------
#
# The variable PYTHON_DIR has to be setup. See ../README.
#

name := exampleB03
G4TARGET := $(name)
G4EXLIB := true


ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk


SWIGINC = include/B03App.i
SWIGCC = src/B03App_wrap.cc
SWIGD = $(G4WORKDIR)/tmp/$(G4SYSTEM)/exampleB03/B03App_wrap.d
objectfilesWrapper = $(G4WORKDIR)/tmp/$(G4SYSTEM)/exampleB03/B03App_wrap.o
SWIGSO = $(G4WORKDIR)/tmp/$(G4SYSTEM)/exampleB03/_B03App.so

moduleName = B03App


MYLIBS = $(subst -l$(name),,$(LDLIBS)) 	

SWIG = $(ANAPHETOP)/specific/$(PLATF)/PublicDomainPackages/2.0.0/bin/swig-1.3.15


.PHONY: all
all: $(SWIGD) $(SWIGCC) lib $(SWIGSO) bin


CPPFLAGS += -I$(ANAPHETOP)/specific/$(PLATF)/PublicDomainPackages/2.0.0/include/python$(PYTHON_VERSION)
CXXFLAGS += -w


LDFLAGS += \
	   -L$(ANAPHETOP)/specific/$(PLATF)/PublicDomainPackages/2.0.0/lib/python$(PYTHON_VERSION)/config \
	   -L/lib \
	   -lnsl -ldl -lreadline /lib/libtermcap.so.2  -lieee -lpthread -lutil -lexampleB03 -lpython$(PYTHON_VERSION)

include $(G4INSTALL)/config/binmake.gmk

$(SWIGSO): $(objectfilesWrapper) 
	@echo Creating shared library:  $@  
	$(CXX) -Wl,-soname,$(@F) -shared $(LDFLAGS) $(objectfilesWrapper) $(MYLIBS) -lc -o $(SWIGSO)

src/%_wrap.cc: include/%.i
	$(SWIG) -python -c++ -shadow -no_default $(SWIGFLAGS) \
	-o $@ $<
	@mv src/$(moduleName).py ./

${G4WORKDIR}/tmp/${G4SYSTEM}/$(name): 
	@mkdir ${G4WORKDIR}/tmp/${G4SYSTEM}/$(name)


#libdir = $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)
#objectfiles = $(shell /bin/ls $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)/*.o)
#_B03App.so:
#	@echo Creating shared library: $(libdir)/$@  
#	$(CXX) -Wl,-soname,$(@F) -shared $(LDFLAGS) $(LDLIBS) -L$(libdir) $(objectfiles) -o $(libdir)/$@ 


