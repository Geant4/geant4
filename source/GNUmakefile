# $Id: GNUmakefile,v 1.1 1999-01-07 16:06:23 gunter Exp $
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

.PHONY: all glob global vis clean_libs clean clean_all

all:
	@for dir in $(SUBDIR1); do (cd $$dir; $(MAKE)); done
	@for dir in $(SUBDIR2); do (cd $$dir; $(MAKE)); done
	@for dir in $(SUBDIR3); do (cd $$dir; $(MAKE)); done

glob global:
	@for dir in $(SUBDIR1); do (cd $$dir; $(MAKE) global); done
	@for dir in $(SUBDIR2); do (cd $$dir; $(MAKE)); done
	@for dir in $(SUBDIR3); do (cd $$dir; $(MAKE)); done

# Following vis target is obsolete and should be removed
#
vis:
	@for dir in $(SUBDIR1); do (cd $$dir; $(MAKE) global); done
	@for dir in $(SUBDIR2); do (cd $$dir; $(MAKE)); done
	@for dir in $(SUBDIR3); do (cd $$dir; $(MAKE)); done

clean_libs:
	@echo Removing all libraries ...
	@for dir in $(SUBDIR1); do (cd $$dir; $(MAKE) clean_libs); done
	@for dir in $(SUBDIR2); do (cd $$dir; $(MAKE) clean_libs); done
	@for dir in $(SUBDIR3); do (cd $$dir; $(MAKE) clean_libs); done
ifdef CFRONT
	@rm -rf $(G4TREP)
endif

clean:
	@echo Removing current $(G4SYSTEM) installation ...
	@rm -rf $(G4TMP)/$(G4SYSTEM)
	@rm -rf $(G4LIB)/$(G4SYSTEM)
	@rm -rf $(G4BIN)/$(G4SYSTEM)
ifdef CFRONT
	@rm -rf $(G4TREP)
endif

# Pay _extremely_ attention to perform the following target !!
#
clean_all:
	@echo Removing all installations ...
	@echo Removing $(G4TMP) ...
	@rm -rf $(G4TMP)
	@echo Removing $(G4LIB) ...
	@rm -rf $(G4LIB)
	@echo Removing $(G4BIN) ...
	@rm -rf $(G4BIN)

