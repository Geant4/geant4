# $Id: GNUmakefile,v 1.80 2006/10/03 13:49:30 gcosmo Exp $
# -----------------------------------------------------------------
# "gmake" makes default libraries for each subdomain.
# "gmake global" makes global libraries for each subdomain.
#                Composite libraries are built.
# "gmake includes" places header files .h/.hh in $G4INCLUDE
# "gmake libmap" forces rebuilding of map-file for granular libraries.
# "gmake dll" forces building of DLLs global libraries and generates DLLs on WIN32
# "gmake clean_libs" removes just archive (.a) and shared (.so) libraries
#                    of current platform installation.
# "gmake clean" removes installation of current platform.
# "gmake clean_all" removes all platforms installations.
#
# (The .o files(s) are made by implicit rules.)

SUBDIR1 =  global intercoms graphics_reps materials
SUBDIR2 =  geometry particles track digits_hits processes
SUBDIR2 += parameterisations tracking event run readout
SUBDIR3 =  persistency interfaces visualization
SUBDIR4 =  g3tog4

ifndef G4INSTALL
  G4INSTALL = ..
endif

include $(G4INSTALL)/config/architecture.gmk
G4BINDIR := $(G4BIN)/$(G4SYSTEM)
G4LIBDIR := $(G4LIB)/$(G4SYSTEM)
G4TMPDIR := $(G4TMP)/$(G4SYSTEM)
unique := $(shell echo $$$$)

.PHONY: all glob global kernel_global kernel libmap dll win32def includes clean_libs clean clean_lists clean_all

all: banner kernel
ifdef G4LIB_USE_G3TOG4
	@for dir in $(SUBDIR4); do (cd $$dir && $(MAKE)); done;:
endif
	@$(MAKE) libmap

glob global: banner kernel_global
ifdef G4LIB_USE_G3TOG4
	@for dir in $(SUBDIR4); do (cd $$dir && $(MAKE)); done;:
endif

kernel_global:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE) global); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) global); done;:
	@for dir in $(SUBDIR3); do (cd $$dir && $(MAKE) global); done;:

kernel:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE)); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE)); done;:
	@for dir in $(SUBDIR3); do (cd $$dir && $(MAKE)); done;:

libmap: $(G4LIBDIR)/liblist
	@$(ECHO) "libmap stage:"
	@$(ECHO) "Searching for GNUmakefiles and sorting ..."
	@find . \
	  -name GNUmakefile -exec $(GREP) -l '^ *name *:=' {} \; \
	  | sort \
	  > /tmp/G4_all_lib_makefiles.$(unique);
	@$(ECHO) "Weeding out paths and files ..."
	@for i in `$(CAT) /tmp/G4_all_lib_makefiles.$(unique)`; \
	do \
	  $(ECHO) $$i | $(GREP) -q -e '/tests\?/' -e 'models/verification' || \
	  $(GREP) -q -e SUBDIR -e G4hepgeometry $$i || \
	  $(ECHO) $$i >> /tmp/G4_granlib_makefiles.$(unique); \
	done
	@$(ECHO) "Making libname.map starter file ..."
	@touch /tmp/G4libname.map.starter.$(unique);
	@for i in `$(CAT) /tmp/G4_granlib_makefiles.$(unique)`; \
	do \
	  $(GREP) '^ *name *:=' $$i | $(CUT) -d \  -f 3 \
	  >> /tmp/G4libname.map.starter.$(unique); \
	  $(ECHO) $$i | $(SED) "s/^\\./source/" \
	  >> /tmp/G4libname.map.starter.$(unique); \
	done
	@$(ECHO) "Making libname.map ..."
	@rm -f $(G4LIBDIR)/libname.map;
	@G4TMP=$(G4TMP); export G4TMP; \
	  $(G4LIBDIR)/liblist -l -d $(G4TMPDIR) \
	  < /tmp/G4libname.map.starter.$(unique) \
	  > $(G4LIBDIR)/libname.map
	@rm -f /tmp/G4_all_lib_makefiles.$(unique);
	@rm -f /tmp/G4_granlib_makefiles.$(unique);
	@rm -f /tmp/G4libname.map.starter.$(unique);

$(G4LIBDIR)/liblist: $(G4INSTALL)/config/liblist.c
	@$(ECHO) "Building library management utility liblist ..."
	@if [ ! -d $(G4LIBDIR) ] ; then mkdir $(G4LIBDIR) ;fi
	@$(CC) $(CCFLAGS) -o $(G4LIBDIR)/liblist $(G4INSTALL)/config/liblist.c

dll: banner win32def
ifneq (,$(findstring WIN32-VC,$(G4SYSTEM)))
	@$(ECHO) "Verifying existence of global libraries for DLLs ..."
	@$(MAKE) kernel_global G4LIB_BUILD_DLL=1
	@$(ECHO) "Building Windows DLL libraries ..."
	@$(MAKE) kernel_global G4LIB_BUILD_DLL=1 G4LIB_BUILD_SHARED=1
	@$(ECHO) "Clearing temporary DLL objects and temporaries ..."
	@rm -rf $(G4LIBDIR)/*.def $(G4LIBDIR)/*.exp
	@rm -rf $(G4TMPDIR)/*
ifdef G4LIB_USE_G3TOG4
	@for dir in $(SUBDIR4); do (cd $$dir && $(MAKE)); done;:
	@mv $(G4LIBDIR)/libG3toG4.a $(G4LIBDIR)/libG3toG4.lib
endif
	@$(ECHO) "Done !"
else
	@$(ECHO) "Sorry !"
endif

win32def: $(G4INSTALL)/config/win32def.c
ifneq (,$(findstring WIN32-VC,$(G4SYSTEM)))
	@$(ECHO) "Building DLL management utility win32def ..."
	@if [ ! -d $(G4LIBDIR) ] ; then mkdir $(G4LIBDIR) ;fi
	@$(CC) $(CCFLAGS) -o $(G4LIBDIR)/win32def $(G4INSTALL)/config/win32def.c
else
	@$(ECHO) "Build of DLLs is only supported on Windows platforms !"
	@$(ECHO) "The system you have currently selected is $(G4SYSTEM)."
endif

includes:
	@$(ECHO) Installing header files in $(G4INCLUDE) ...
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE) $@); done
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) $@); done
	@for dir in $(SUBDIR3); do (cd $$dir && $(MAKE) $@); done
ifdef G4LIB_USE_G3TOG4
	@for dir in $(SUBDIR4); do (cd $$dir && $(MAKE) $@); done
endif

banner:
	@$(ECHO) "*************************************************************"
	@$(ECHO) " Installation Geant4 version $$Name: geant4-08-01-patch-02 $ "
	@$(ECHO) " Copyright (C) 1994-2006 Geant4 Collaboration                            "
	@$(ECHO) "*************************************************************"

clean_libs::
	@$(ECHO) Removing all libraries ...
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE) clean_libs); done
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) clean_libs); done
	@for dir in $(SUBDIR3); do (cd $$dir && $(MAKE) clean_libs); done
ifdef G4LIB_USE_G3TOG4
	@for dir in $(SUBDIR4); do (cd $$dir && $(MAKE) clean_libs); done
endif
	@rm -rf $(G4LIBDIR)/libname.map

clean::
	@$(ECHO) Removing current $(G4SYSTEM) installation ...
	@rm -rf $(G4TMPDIR)
	@rm -rf $(G4LIBDIR)
	@rm -rf $(G4BINDIR)
	@rm -rf $(G4INSTALL)/.config/bin/$(G4SYSTEM)
	@$(ECHO) Removing physics-lists installation ...
ifndef G4LISTS_BASE
	@rm -rf $(G4LIB)/plists/$(G4SYSTEM)
	@rm -rf $(G4TMP)/plists/$(G4SYSTEM)
else
	@rm -rf $(G4LISTS_BASE)/hadronic/plists/lib/$(G4SYSTEM)
	@rm -rf $(G4LISTS_BASE)/hadronic/plists/tmp/$(G4SYSTEM)
endif
ifdef CFRONT
	@rm -rf $(G4TREP)
endif

clean_lists::
	@$(ECHO) Removing physics-lists installation ...
ifndef G4LISTS_BASE
	@rm -rf $(G4LIB)/plists/$(G4SYSTEM)
	@rm -rf $(G4TMP)/plists/$(G4SYSTEM)
else
	@rm -rf $(G4LISTS_BASE)/hadronic/plists/lib/$(G4SYSTEM)
	@rm -rf $(G4LISTS_BASE)/hadronic/plists/tmp/$(G4SYSTEM)
endif

# Pay _extremely_ attention before executing the following target !!
#
clean_all::
	@$(ECHO) Removing all installations ...
	@$(ECHO) Removing $(G4TMP) ...
	@rm -rf $(G4TMP)
	@$(ECHO) Removing $(G4LIB) ...
	@rm -rf $(G4LIB)
	@$(ECHO) Removing $(G4BIN) ...
	@rm -rf $(G4BIN)
	@$(ECHO) Removing $(G4INCLUDE) ...
	@rm -rf $(G4INCLUDE)
ifdef G4LISTS_BASE
	@$(ECHO) Removing physics-lists installation ...
	@rm -rf $(G4LISTS_BASE)/hadronic/plists
endif
	@$(ECHO) Removing installation configurations ...
	@rm -rf $(G4INSTALL)/.config
