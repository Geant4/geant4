# $Id: GNUmakefile,v 1.12 1999-12-14 13:58:46 stesting Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

ifndef G4INSTALL
  G4INSTALL = ..
endif
 
include $(G4INSTALL)/config/architecture.gmk

#Genuine tests...
SUBDIRS  = test01 ## Geometry - CSG solids and navigator.
                  ## (Action? Augment, plus Booleans.)
SUBDIRS += test02 ## Standard em processes.  Low energy and high energy
                  ## hadronic processes and stopping, including coherent
                  ## elastic scattering.  (Derived originally from ExN04.)
                  ## (Action? Split into (a) em, (b) had with higher em cuts.)
# SUBDIRS += test03 ## Optical Processes.  (Derived originally from ExN06.)
                    ## Archived (cvs rm) - was indentical to ExN06 and no
                    ## longer maintained.
# SUBDIRS += test04 ## Persistency.
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
SUBDIRS += test17 ## LowEnergy em (protons, antiprotons and ions).  (Derived
                  ## from TestEm6.)

# Tests with hadronic exercisers.
HADEXERDIRS = test02 test11 test12 test13 test15 test16

# Symbolic links to novice examples...
EXAMPLES = 01 02 03 04 05 06

# Symbolic links to persistency examples...
PERSISTENCY += 01 02

.PHONY : all examples_links clean clean_libs clean_examples_links

all:
	@$(ECHO) Making libG4TestTool.a...
	@cd tools/lib; $(MAKE)
	@$(ECHO) Making hadronic exercisers...
	@for dir in $(HADEXERDIRS); \
	  do (cd $$dir; $(MAKE) -f GNUmakefile.hadronic.exerciser); done
	@$(ECHO) Making tests...
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done
	@$(ECHO) Making links to examples...
	@$(MAKE) examples_links
	@$(ECHO) Making examples...
	@for nn in $(EXAMPLES); do \
	  (cd test1$$nn; $(MAKE)); \
	done
ifdef G4USE_HEPODBMS
	@for nn in $(PERSISTENCY); do \
	  (cd test4$$nn; $(MAKE)); \
	done
endif

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done
	@for nn in $(EXAMPLES); do \
	  (cd test1$$nndir; $(MAKE) clean); \
	done
	@$(MAKE) clean_examples_links

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
	@for nn in $(EXAMPLES); do \
	  (cd test1$$nn; $(MAKE) clean_libs); \
	done
	@$(MAKE) clean_examples_links

examples_links:
	@$(MAKE) clean_examples_links
	@for nn in $(EXAMPLES); do \
	  testdir=test1$$nn; \
	  testin=test1$$nn.in; \
	  testout=test1$$nn.out; \
	  exdir=$(G4INSTALL)/examples/novice/N$$nn; \
	  exin=exampleN$$nn.in; \
	  exout=exampleN$$nn.out; \
	  ln -s $$exdir $$testdir; \
	  ln -s $$exdir/$$exin $$testdir/$$testin; \
	  ln -s $$exdir/$$exout $$testdir/$$testout; \
	  ln -s $(G4WORKDIR)/bin/$(G4SYSTEM)/exampleN$$nn \
	    $(G4WORKDIR)/bin/$(G4SYSTEM)/test1$$nn; \
	done
	@ln -s $(G4INSTALL)/examples/novice/N04/exampleN04.EMtest.in \
          test104/test104.EMtest.in
	@ln -s $(G4INSTALL)/examples/novice/N04/exampleN04.EMtest.out \
          test104/test104.EMtest.out
	@for nn in $(PERSISTENCY); do \
	  testdir=test4$$nn; \
	  testin=test4$$nn.in; \
	  testout=test4$$nn.out; \
	  exdir=$(G4INSTALL)/examples/extended/persistency/PersistentEx$$nn; \
	  exin=PersistentEx$$nn.in; \
	  exout=PersistentEx$$nn.out; \
	  ln -s $$exdir $$testdir; \
	  ln -s $$exdir/$$exin $$testdir/$$testin; \
	  ln -s $$exdir/$$exout $$testdir/$$testout; \
	  ln -s $$exdir/PersistentEx$${nn}run \
	    $(G4WORKDIR)/bin/$(G4SYSTEM)/test4$$nn; \
	done

clean_examples_links:
	@rm -f test104/test104.EMtest.in
	@rm -f test104/test104.EMtest.out
	@for nn in $(EXAMPLES); do \
	  testdir=test1$$nn; \
	  testin=test1$$nn.in; \
	  testout=test1$$nn.out; \
	  exdir=$(G4INSTALL)/examples/novice/N$$nn; \
	  rm -f $$testdir/$$testin; \
	  rm -f $$testdir/$$testout; \
	  rm -f $$testdir; \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test1$$nn; \
	done
	@for nn in $(PERSISTENCY); do \
	  testdir=test4$$nn; \
	  testin=test4$$nn.in; \
	  testout=test4$$nn.out; \
	  exdir=$(G4INSTALL)/examples/extended/persistency/PersistentEx$$nn; \
	  rm -f $$testdir/$$testin; \
	  rm -f $$testdir/$$testout; \
	  rm -f $$testdir; \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test4$$nn; \
	done
