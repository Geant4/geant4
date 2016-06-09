# $Id: GNUmakefile,v 1.12 2004/06/20 17:26:13 daquinog Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------
#
# The variable PI_BASE_DIR has to be setup. See ../README.
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

SWIG = ${SWIG_BASE_DIR}/bin/swig-${SWIGVERSION}

.PHONY: clean_all all
all: $(SWIGD) $(SWIGCC) lib $(SWIGSO) bin


CPPFLAGS += -I${PI_BASE_DIR}/${PI_VER}/external/Python/2.2.2/${PLATF}/include/python${PYTHONVERSION}
CXXFLAGS += -w

LDFLAGS += \
	   -L${PI_BASE_DIR}/${PI_VER}/external/Python/2.2.2/${PLATF}/lib/python${PYTHONVERSION}/config \
	   -L/lib \
	   -lnsl -ldl -lreadline /lib/libtermcap.so.2  -lieee -lpthread -lutil -lexampleB03 -lpython$(PYTHONVERSION)

include $(G4INSTALL)/config/binmake.gmk

ifdef G4ANALYSIS_USE
 CPPFLAGS += `aida-config --include`
 LDFLAGS  += `aida-config --lib`
endif

$(SWIGSO): $(objectfilesWrapper) 
	@echo Creating shared library:  $@  
	$(CXX) -Wl,-soname,$(@F) -shared $(LDFLAGS) $(objectfilesWrapper) $(MYLIBS) -lc -o $(SWIGSO)

src/%_wrap.cc: include/%.i
	$(SWIG) -python -c++ -shadow -no_default $(SWIGFLAGS) \
	-o $@ $<
	@mv src/$(moduleName).py ./

${G4WORKDIR}/tmp/${G4SYSTEM}/$(name): 
	@mkdir ${G4WORKDIR}/tmp/${G4SYSTEM}/$(name)


clean_all: clean
	@rm -f B03App.pyc
	@rm B03App.py
	@rm src/B03App_wrap.cc


#libdir = $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)
#objectfiles = $(shell /bin/ls $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(name)/*.o)
#_B03App.so:
#	@echo Creating shared library: $(libdir)/$@  
#	$(CXX) -Wl,-soname,$(@F) -shared $(LDFLAGS) $(LDLIBS) -L$(libdir) $(objectfiles) -o $(libdir)/$@ 


