# $Id: GNUmakefile,v 1.26 2000-06-15 16:21:26 stesting Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

ifndef G4INSTALL
  G4INSTALL = ..
endif
 
include $(G4INSTALL)/config/architecture.gmk

# Symbolic links to novice examples using test101, etc...
EXAMPLES = 01 02 03 04 05 06

# Symbolic links to extended electromagnetic examples using test501, etc...
ELECTROMAGNETIC += 1 2 3 4 5 ## Awaiting description.
# ELECTROMAGNETIC += 6 ## Low energy electromagnetic processes (moved to lowenergy).
# ELECTROMAGNETIC += 8 ## Transistion radiation and PAI (photoabsorption
                       ## ionisation) model.

# Symbolic links to extended persistency examples using test401, etc...
PERSISTENCY += 01 02

# System Tests...
SUBDIRS += test01 ## Geometry - CSG solids and navigator.
                  ## (Action? Augment, plus Booleans.)
SUBDIRS += test02 ## Standard em processes.  Low energy and high energy
                  ## hadronic processes and stopping, including coherent
                  ## elastic scattering.  (Derived originally from ExN04.)
                  ## (Action? Split into (a) em, (b) had with higher em cuts.)
# SUBDIRS += test03 ## Optical Processes.  (Derived originally from ExN06.)
                    ## Archived (cvs rm) - was indentical to ExN06 and no
                    ## longer maintained.
# SUBDIRS += test04 ## Persistency.  Retired.
SUBDIRS += test05 ## BREPS and STEP Interface.
SUBDIRS += test06 ## Integral approach electromagnetic processes.  Simple box.
                  ## (Action?  Add muon test script.)
SUBDIRS += test07 ## Statistical evaluation prototype.  (Derived originally
                  ## from ExN03.)
# SUBDIRS += test08 ## Old em test.  Archived (cvs rm) for now.
                    ## (Action? Replace with transition radiation
                    ## and PAI (photoabsorption ionisation) model.)  (Derived
                    ## from TestEm8.)
SUBDIRS += test09 ## Basic e.m. processes, Polarised Compton, new photelectric
                  ## effect, synchrotron radiation, simple box.  (Action?  Add
                  ## magnetic field, ElossPLus, muon test script.)
SUBDIRS += test10 ## Uses optical photons to test geometry.
SUBDIRS += test11 ## Neutron transport
SUBDIRS += test12 ## FTF String Model
SUBDIRS += test13 ## Dual Parton Model
SUBDIRS += test14 ## LowEnergy em (photons and electrons).  (Derived from
                  ## TestEm5.)
SUBDIRS += test15 ## Isotope Production
SUBDIRS += test16 ## n and p Cross-Sections
SUBDIRS += test17 ## LowEnergy em (protons, antiprotons and ions). 
SUBDIRS += test18 ## Radioactive decay.

# Tests with hadronic exercisers.
HADEXERDIRS = test02 test11 test12 test13 test15 test16

.PHONY : all examples_links clean clean_libs clean_examples_links

all:
# Making links in the source tree can be a problem if more than one
# jobs try to do it at the same time.  So disable and always do by
# hand with: gmake examples_links_source
	@$(ECHO) Making examples...
	@for nn in $(EXAMPLES); do \
	  (cd test1$$nn && $(MAKE)); \
	done;:
	@for nn in $(ELECTROMAGNETIC); do \
	  (cd test50$$nn && $(MAKE)); \
	done;:
ifdef G4USE_HEPODBMS
	@for nn in $(PERSISTENCY); do \
	  (cd test4$$nn && $(MAKE)); \
	done;:
endif
	@$(ECHO) Making libG4TestTool.a for System Tests...
	@cd tools/lib && $(MAKE);:
	@$(ECHO) Making hadronic exercisers for System Tests...
	@for dir in $(HADEXERDIRS); \
	  do (cd $$dir && $(MAKE) -f GNUmakefile.hadronic.exerciser); done;:
	@$(ECHO) Making System Tests...
	@for dir in $(SUBDIRS); do (cd $$dir && $(MAKE)); done;:
	@$(ECHO) Making links to examples binaries...
	@$(MAKE) examples_links_binaries

clean:
	@for nn in $(EXAMPLES); do \
	  (cd test1$$nndir && $(MAKE) clean); \
	done;:
	@for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) clean); done;:

clean_libs:
	@for nn in $(EXAMPLES); do \
	  (cd test1$$nn && $(MAKE) clean_libs); \
	done;:
	@for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) clean_libs); done;:

examples_links_source:
	@$(MAKE) clean_examples_links_source
	@for nn in $(EXAMPLES); do \
	  testdir=test1$$nn; \
	  exdir=$(G4INSTALL)/examples/novice/N$$nn; \
	  ln -s $$exdir $$testdir; \
	  for i in $$exdir/*.in; do \
	    exin=`basename $$i`; \
	    exout=`basename $$i .in`.out; \
	    testin=`echo $$exin | sed s/exampleN/test1/`; \
	    testout=`echo $$exout | sed s/exampleN/test1/`; \
	    ln -s $$exin $$testdir/$$testin; \
	    ln -s $$exout $$testdir/$$testout; \
	  done; \
	done;:
	@for nn in $(ELECTROMAGNETIC); do \
	  testdir=test50$$nn; \
	  exdir=$(G4INSTALL)/examples/extended/electromagnetic/TestEm$$nn; \
	  ln -s $$exdir $$testdir; \
	  for i in $$exdir/*.in; do \
	    exin=`basename $$i`; \
	    exout=`basename $$i .in`.out; \
	    testin=`echo $$exin | sed s/TestEm/test50/`; \
	    testout=`echo $$exout | sed s/TestEm/test50/`; \
	    ln -s $$exin $$testdir/$$testin; \
	    ln -s $$exout $$testdir/$$testout; \
	  done; \
	done;:
	@for nn in $(PERSISTENCY); do \
	  testdir=test4$$nn; \
	  exdir=$(G4INSTALL)/examples/extended/persistency/PersistentEx$$nn; \
	  ln -s $$exdir $$testdir; \
	  for i in $$exdir/*.in; do \
	    exin=`basename $$i`; \
	    exout=`basename $$i .in`.out; \
	    testin=`echo $$exin | sed s/PersistentEx/test4/`; \
	    testout=`echo $$exout | sed s/PersistentEx/test4/`; \
	    ln -s $$exin $$testdir/$$testin; \
	    ln -s $$exout $$testdir/$$testout; \
	  done; \
	done;:

examples_links_binaries:
	@$(MAKE) clean_examples_links_binaries
	@for nn in $(EXAMPLES); do \
	  ln -s exampleN$$nn $(G4WORKDIR)/bin/$(G4SYSTEM)/test1$$nn; \
	done;:
	@for nn in $(ELECTROMAGNETIC); do \
	  ln -s TestEm$$nn $(G4WORKDIR)/bin/$(G4SYSTEM)/test50$$nn; \
	done;:
	@for nn in $(PERSISTENCY); do \
	  exdir=$(G4INSTALL)/examples/extended/persistency/PersistentEx$$nn; \
	  ln -s $$exdir/PersistentEx$${nn}run \
	    $(G4WORKDIR)/bin/$(G4SYSTEM)/test4$$nn; \
	done;:

clean_examples_links_source:
	@rm -f test104/test104.EMtest.in
	@rm -f test104/test104.EMtest.out
	@for nn in $(EXAMPLES); do \
	  testdir=test1$$nn; \
	  for i in $$testdir/$$testdir.*; do \
	    rm -f $$i; \
	  done; \
	  rm -f $$testdir; \
	done;:
	@for nn in $(ELECTROMAGNETIC); do \
	  testdir=test50$$nn; \
	  for i in $$testdir/$$testdir.*; do \
	    rm -f $$i; \
	  done; \
	  rm -f $$testdir; \
	done;:
	@for nn in $(PERSISTENCY); do \
	  testdir=test4$$nn; \
	  for i in $$testdir/$$testdir.*; do \
	    rm -f $$i; \
	  done; \
	  rm -f $$testdir; \
	done;:

clean_examples_links_binaries:
	@for nn in $(EXAMPLES); do \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test1$$nn; \
	done;:
	@for nn in $(ELECTROMAGNETIC); do \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test50$$nn; \
	done;:
	@for nn in $(PERSISTENCY); do \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test4$$nn; \
	done;:
