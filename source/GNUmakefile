# $Id: GNUmakefile 104933 2017-06-30 08:07:44Z gcosmo $
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

ifndef G4INSTALL
  G4INSTALL = ..
endif

include $(G4INSTALL)/config/architecture.gmk

SUBDIR1 =  externals global analysis intercoms graphics_reps materials
SUBDIR2 =  geometry particles track digits_hits processes parameterisations
SUBDIR2 += tracking event run error_propagation readout physics_lists
SUBDIR3 =  persistency interfaces visualization

SUBDIR4 =  g3tog4

G4BINDIR := $(G4BIN)/$(G4SYSTEM)
G4LIBDIR := $(G4LIB)/$(G4SYSTEM)
G4TMPDIR := $(G4TMP)/$(G4SYSTEM)
unique := $(shell echo $$$$)

.PHONY: all glob global obj kernel_global kernel liblist libmap dll win32def includes clean_libs clean clean_lists clean_all

all: banner kernel
ifdef G4LIB_USE_G3TOG4
	@for dir in $(SUBDIR4); do (cd $$dir && $(MAKE)); done;:
endif
	@$(MAKE) libmap
	@$(ECHO) "Libraries installation completed !"

glob global: banner kernel_global
ifdef G4LIB_USE_G3TOG4
	@for dir in $(SUBDIR4); do (cd $$dir && $(MAKE)); done;:
endif
	@$(MAKE) liblist
	@$(ECHO) "Libraries installation completed !"

fpe_debug: banner kernel
	@$(ECHO) Activating FPE debug build ...
ifneq (,$(findstring Linux-g++,$(G4SYSTEM)))
	@rm -f $(G4TMPDIR)/G4run/G4RunManagerKernel.*
	@cd $(G4INSTALL)/source/run; $(MAKE) G4FPE_DEBUG=1;
	@cd $(G4INSTALL)/source; $(MAKE) libmap;
	@$(ECHO) "Libraries installation completed !"
else
	@$(ECHO) "Build of FPE check is only supported on Linux platforms with gcc !"
	@$(ECHO) "The system you have currently selected is $(G4SYSTEM)."
	@$(ECHO) "No actions taken !"
endif

kernel_global:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE) global); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) global); done;:
	@for dir in $(SUBDIR3); do (cd $$dir && $(MAKE) global); done;:

kernel:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE)); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE)); done;:
	@for dir in $(SUBDIR3); do (cd $$dir && $(MAKE)); done;:

obj: banner
	@$(ECHO) Building all object files first ...
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE) obj); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) obj); done;:
	@for dir in $(SUBDIR3); do (cd $$dir && $(MAKE) obj); done;:
ifdef G4LIB_USE_G3TOG4
	@for dir in $(SUBDIR4); do (cd $$dir && $(MAKE) obj); done;:
endif

libmap: liblist
	@$(ECHO) "Libmap stage. Searching for GNUmakefiles and sorting ..."
	@$(FIND) . \
	  -name GNUmakefile -exec $(GREP) -l '^ *name *:=' {} \; \
	  | $(SORT) \
	  > /tmp/G4_all_lib_makefiles.$(unique);
	@$(ECHO) "Weeding out paths and files ..."
	@for i in `$(CAT) /tmp/G4_all_lib_makefiles.$(unique)`; \
	do \
	  $(ECHO) $$i | $(GREP) -q -e '/tests\?/' -e 'models/verification' || \
	  $(GREP) -q -e SUBDIR -e G4hepgeometry $$i || \
	  $(ECHO) $$i >> /tmp/G4_granlib_makefiles.$(unique); \
	done
	@$(ECHO) "Making libname.map starter file ..."
	@$(TOUCH) /tmp/G4libname.map.starter.$(unique);
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

liblist: $(G4LIBDIR)/liblist

$(G4LIBDIR)/liblist: $(G4INSTALL)/config/liblist.c
	@$(ECHO) "Building library management utility liblist ..."
	@if [ ! -d $(G4LIBDIR) ] ; then mkdir $(G4LIBDIR) ;fi
	@$(CC) $(CCFLAGS) -o $(G4LIBDIR)/liblist $(G4INSTALL)/config/liblist.c

clean_tmps:
	@$(ECHO) "Clearing temporary objects ..."
	@rm -rf $(G4LIBDIR)/*.def $(G4LIBDIR)/*.exp
	@rm -rf $(G4TMPDIR)/*

dll_build: banner genwindef
ifneq (,$(findstring WIN32-VC,$(G4SYSTEM)))
	@$(ECHO) "Verifying existence of global libraries for DLLs ..."
	@$(MAKE) kernel_global G4LIB_BUILD_DLL=1
ifdef G4LIB_USE_G3TOG4
	@for dir in $(SUBDIR4); do (cd $$dir && $(MAKE) G4LIB_BUILD_DLL=1); done;:
endif
	@$(ECHO) "Building Windows DLL libraries ..."
	@$(MAKE) kernel_global G4LIB_BUILD_DLL=1 G4LIB_BUILD_SHARED=1
ifdef G4LIB_USE_G3TOG4
	@for dir in $(SUBDIR4); do (cd $$dir && $(MAKE) G4LIB_BUILD_DLL=1 G4LIB_BUILD_SHARED=1); done;:
endif
	@$(ECHO) "Done !"
endif

dll:
ifneq (,$(findstring WIN32-VC,$(G4SYSTEM)))
	@$(MAKE) dll_build
	@$(MAKE) clean_tmps
else
	@$(ECHO) "Sorry !"
endif

genwindef: $(G4INSTALL)/config/genwindef.cc
ifneq (,$(findstring WIN32-VC,$(G4SYSTEM)))
	@$(ECHO) "Building DLL management utility genwindef ..."
	@if [ ! -d $(G4LIBDIR) ] ; then mkdir $(G4LIBDIR) ;fi
	@$(CXX) $(CXXFLAGS) /Fe$(G4LIBDIR)/genwindef.exe $(G4INSTALL)/config/genwindef.cc
else
	@$(ECHO) "Build of DLLs is only supported on Windows platforms !"
	@$(ECHO) "The system you have currently selected is $(G4SYSTEM)."
	@$(ECHO) "No actions taken !"
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
	@$(ECHO) " Installation Geant4 version geant4-10-04-beta-01            "
	@$(ECHO) " Copyright (C) 1994-2017 Geant4 Collaboration                "
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
	@rm -rf $(G4INSTALL)/source/interfaces/*/moc
	@rm -rf $(G4INSTALL)/source/visualization/*/moc
	@rm -rf $(G4INSTALL)/.config/bin/$(G4SYSTEM)
ifdef CFRONT
	@rm -rf $(G4TREP)
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
	@$(ECHO) Removing Qt moc source files ...
	@rm -rf $(G4INSTALL)/source/interfaces/*/moc
	@rm -rf $(G4INSTALL)/source/visualization/*/moc
	@$(ECHO) Removing installation configurations ...
	@rm -rf $(G4INSTALL)/.config
