# $Id:$
# ------------------------------------------------------------
# GNUmakefile for internal memoryprotection module.
# ------------------------------------------------------------

G4TARGET := tracer

name := G4memoryprotection

ifndef G4INSTALL
  G4INSTALL := ../..
endif

GLOBLIBS  = libG4global.lib

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS  += -I$(G4BASE)/global/management/include

.SUFFIXES: 
.SUFFIXES:	.a .o .cc .cxx .c .d

G4LIBDIR := $(G4LIB)/$(G4SYSTEM)
G4TMPDIR := $(G4TMP)/$(G4SYSTEM)/$(name)
G4BINDIR := $(G4BIN)/$(G4SYSTEM)

Tbinsources = src/$(G4TARGET).cxx
Tbinobjects = $(patsubst src/%.cxx,$(G4TMPDIR)/%.o,$(Tbinsources))
Tdependencies = $(patsubst %.o,%.d,$(Tbinobjects))

.PHONY: global all lib bin libc

###--------------------- section to define explicit rules --------------------

lib: libc

global: cleanf lib bin

all: lib

ifneq ($(G4LIB_BUILD_SHARED),)
libc : $(G4LIBDIR)/lib$(name).$(SHEXT)
else
libc : $(G4LIBDIR)/lib$(name).a
endif

bin: $(G4BINDIR)/$(G4TARGET)
	@echo $@ stage done

### Standalone program for tracing memory
#
$(G4BINDIR)/$(G4TARGET) : $(Tbinobjects)
	@if [ ! -d $(G4BINDIR) ] ; then mkdir $(G4BINDIR) ;fi
	@echo Linking $(G4TARGET) ...
	@$(RM) $@
ifneq ($(G4LIB_BUILD_SHARED),)
	@$(CC) $(Tbinobjects) -o $@ -lpthread -lstdc++ $(G4LIBDIR)/lib$(name).$(SHEXT)
else
	@$(CC) $(Tbinobjects) -o $@ -lpthread -lstdc++ $(G4LIBDIR)/lib$(name).a
endif
	@echo Finished linking $(G4TARGET)

###-------------------- locally defined pattern rules ------------------------

### Make the .o files from src/*.cxx
#
$(G4TMPDIR)/%.o : src/%.cxx libc
	@echo Building $(G4TARGET) ...
	@$(CC) $(CXXFLAGS) $(CPPFLAGS) -o $(G4TMPDIR)/$(G4TARGET).o -c src/$(G4TARGET).cxx

### --------------------- .d files --------------------------------------------

### Make the .d file(s) and include it(them).
#
$(G4TMPDIR)/%.d: src/%.cxx
	@if [ ! -d $(G4TMP)/$(G4SYSTEM) ] ; then mkdir $(G4TMP)/$(G4SYSTEM) ;fi
	@if [ ! -d $(G4TMPDIR) ] ; then mkdir $(G4TMPDIR) ;fi
	@if test -f $@ ; then rm -f $@; fi

### --------------------- end of pattern rules --------------------------------

clean:: cleanf

cleanf:
	@$(RM) $(G4BINDIR)/$(G4TARGET)

include $(G4INSTALL)/config/common.gmk
