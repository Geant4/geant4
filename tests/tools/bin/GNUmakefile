include $(G4INSTALL)/config/architecture.gmk

DUMMY := $(shell if [ ! -d $(G4SYSTEM) ]\; then mkdir $(G4SYSTEM) \; fi)

.PHONY: $(G4SYSTEM)/.depend

vpath %.hh ../include
vpath %.cc ../binsrc

CPPFLAGS+=-I../include

AR:=ar

%: %.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS)  $< $(LDFLAGS) -o$@


SRCS:=$(notdir $(wildcard ../binsrc/*.cc))
BINS:=$(patsubst %.cc,%,$(SRCS))
ifeq ($(G4SYSTEM),G++-Linux)
EXTRADEPEND=-I/usr/include/g++
endif

LDFLAGS+=../lib/$(G4SYSTEM)/libG4TestTool.a

all: $(BINS)
	mv $(BINS) $(G4SYSTEM)

$(G4SYSTEM)/.depend:
	if [ ! -f $(G4SYSTEM)/.depend ]; then touch $(G4SYSTEM)/.depend ; fi
	cd ../binsrc; \
        makedepend -f ../bin/$(G4SYSTEM)/.depend -I../include $(EXTRADEPEND) $(SRCS)	

include $(G4SYSTEM)/.depend





