# $Id: GNUmakefile,v 1.7 1999-05-07 13:13:40 stesting Exp $
# -----------------------------------------------------------------
# "gmake" makes default libraries for each subdomain.
# "gmake glob" makes global libraries for each subdomain.
#         Composite libraries are built.
# "gmake vis" makes global libraries and default (separate) visualization.
# "gmake total" makes everything, including the BREPS and the STEP code.
#
# (The .o files(s) are made by implicit rules.)

MAKEFLAGS= --no-print-directory

SUBDIR1 = global particles processes geometry digits+hits
SUBDIR2 = event tracking track materials run graphics_reps intercoms readout
SUBDIR3 = g3tog4 interfaces visualization
SUBDIR4 = # persistency

ifndef G4INSTALL
  G4INSTALL = ..
endif

include $(G4INSTALL)/config/architecture.gmk
G4BINDIR := $(G4BIN)/$(G4SYSTEM)
unique := $(shell echo $$$$)

.PHONY: all glob global libmap clean_libs clean clean_all

all:
	@for dir in $(SUBDIR1); do (cd $$dir; $(MAKE)); done
	@for dir in $(SUBDIR2); do (cd $$dir; $(MAKE)); done
	@for dir in $(SUBDIR3); do (cd $$dir; $(MAKE)); done
	$(MAKE) libmap

glob global:
	@for dir in $(SUBDIR1); do (cd $$dir; $(MAKE) global); done
	@for dir in $(SUBDIR2); do (cd $$dir; $(MAKE)); done
	@for dir in $(SUBDIR3); do (cd $$dir; $(MAKE)); done

libmap: $(G4BINDIR)/liblist
	@echo "WARNING: Making a library map of granular libraries."
	@echo "         This is a list of libraries in order of use, and for"
	@echo "         each library a list of other libraries used."
	@echo "         To do this is needs a complete set of dependency"
	@echo "         files, e.g., after gmake in the source/ directory."
	@echo "Searching $(G4INSTALL)/source"
	@echo '  for GNUmakefiles containing "name" and sorting...'
	@find $(G4INSTALL)/source \
	  -name GNUmakefile -exec $(GREP) -l 'name :=' {} \; \
	  | sort \
	  > /tmp/G4_all_lib_makefiles.$(unique);
	@echo "Weeding out global level GNUmakefiles and non-libraries..."
	@for i in `cat /tmp/G4_all_lib_makefiles.$(unique)`; \
	do \
	  $(GREP) -q SUBDIR $$i || \
	  $(GREP) -q G4hepgeometry $$i || \
	  echo $$i >> /tmp/G4_granlib_makefiles.$(unique); \
	done
	@echo "Making libname.map starter file..."
	@touch /tmp/G4libname.map.starter.$(unique);
	@for i in `cat /tmp/G4_granlib_makefiles.$(unique)`; \
	do \
	  $(GREP) 'name :=' $$i | cut -d \  -f 3 \
	  >> /tmp/G4libname.map.starter.$(unique); \
	  echo $$i >> /tmp/G4libname.map.starter.$(unique); \
	done
	@echo "Making libname.map..."
	@if [ ! -d $(G4TMP)/$(G4SYSTEM) ]; then mkdir $(G4TMP)/$(G4SYSTEM);fi
	@rm -f $(G4TMP)/$(G4SYSTEM)/libname.map;
	$(G4BINDIR)/liblist -l \
	  -d $(G4TMP)/$(G4SYSTEM) \
	  < /tmp/G4libname.map.starter.$(unique) \
	  > $(G4TMP)/$(G4SYSTEM)/libname.map
	@rm -f /tmp/G4_all_lib_makefiles.$(unique);
	@rm -f /tmp/G4_granlib_makefiles.$(unique);
	@rm -f /tmp/G4libname.map.starter.$(unique);

$(G4BINDIR)/liblist: $(G4INSTALL)/config/liblist.c
	@echo "Compiling liblist.c..."
	@if [ ! -d $(G4BINDIR) ] ; then mkdir $(G4BINDIR) ;fi
	$(CC) $(CCFLAGS) -o $(G4BINDIR)/liblist $(G4INSTALL)/config/liblist.c

clean_libs:
	@echo Removing all libraries ...
	@for dir in $(SUBDIR1); do (cd $$dir; $(MAKE) clean_libs); done
	@for dir in $(SUBDIR2); do (cd $$dir; $(MAKE) clean_libs); done
	@for dir in $(SUBDIR3); do (cd $$dir; $(MAKE) clean_libs); done
ifdef CFRONT
	@rm -rf $(G4TREP)
endif
	@rm -rf $(G4TMP)/$(G4SYSTEM)/libname.map

clean:
	@echo Removing current $(G4SYSTEM) installation ...
	@rm -rf $(G4TMP)/$(G4SYSTEM)
	@rm -rf $(G4LIB)/$(G4SYSTEM)
	@rm -rf $(G4BIN)/$(G4SYSTEM)
ifdef CFRONT
	@rm -rf $(G4TREP)
endif

# Pay _extremely_ attention before executing the following target !!
#
clean_all:
	@echo Removing all installations ...
	@echo Removing $(G4TMP) ...
	@rm -rf $(G4TMP)
	@echo Removing $(G4LIB) ...
	@rm -rf $(G4LIB)
	@echo Removing $(G4BIN) ...
	@rm -rf $(G4BIN)
