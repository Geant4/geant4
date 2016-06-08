include $(HEP_ODBMS_DIR)/etc/HepODBMS.mk

# application name, source files on objects
APPL       = createTag
SRCS       = $(APPL).cpp
OBJS	   = $(OS)/$(APPL).$(OBJ_EXT)
LIBS       = $(HEP_ODBMS_LIBS) $(OBJY_LIBS) $(OSPACE_LIBS)

G4INSTALL         := ../../../../..
G4ODBMS           := true
ifeq ($(OS),Solaris)
  G4SYSTEM        := SUN-CC
endif
ifeq ($(OS),HP-UX)
  G4SYSTEM        := HP-aCC
endif

CLHEPLIB   := /afs/cern.ch/sw/geant4/dev/CLHEP/$(G4SYSTEM)/pro/lib/libCLHEP.a

include $(G4INSTALL)/config/architecture.gmk

.PHONY: cleanall

# CPPFLAGS is defined in Geant4 architecture.gmk
CPPFLAGS  += \
	     -I$./ \
             -I$(G4SCHEMA_INCLUDE) \
	     -I$(G4INSTALL)/source/global/management/include \
             -I$(G4INSTALL)/source/persistency/global/include \
             -I$(G4INSTALL)/source/persistency/management/include \
             -I$(G4INSTALL)/source/persistency/events/include \
             -I$(G4INSTALL)/source/persistency/digits+hits/hits/include \
             -I$(G4INSTALL)/source/event/include \
             -I$(G4INSTALL)/source/tracking/include \
             -I$(G4INSTALL)/source/digits+hits/hits/include \
             -I$(G4INSTALL)/source/digits+hits/digits/include

# C_FLAGS is the one used in the rules of HepODBMS.mk
C_FLAGS += $(CPPFLAGS)

# Add Geant4 Persistent Libraries
LIBS :=      -L$(G4LIB)/$(G4SYSTEM) \
             -lG4pmanagement -lG4pgeomn -lG4pgeomBoolean  -lG4pcsg \
             -lG4pgeomGlobal -lG4prun -lG4pevents -lG4phits -lG4pglobal \
             -lG4event -lG4partman -lG4hits -lG4digits \
             -lG4materials \
             -lG4intercoms -lG4globman -lG4hepnumerics \
             $(LIBS) \
             $(CLHEPLIB)

# make APPL the default target
all : $(APPL)

$(APPL): GNUmakefile $(OBJS)
	$(C++) $(C_FLAGS) $(CPPFLAGS) -o $(APPL) $(OBJS) $(LIBS)

cleanall : clean
	@rm -f $(APPL)

