# $Id: GNUmakefile,v 1.11 1999-12-10 11:50:33 stesting Exp $
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

# Symbolic links to examples...
EXAMPLES = 01 02 03 04 05 06

.PHONY : all clean clean_libs

all:
	@rm -f test104/test104.EMtest.in
	@rm -f test104/test104.EMtest.out
	@rm -f test104
	@ln -s $(G4INSTALL)/examples/novice/N04 test104
	@ln -s $(G4INSTALL)/examples/novice/N04/exampleN04.EMtest.in \
          test104/test104.EMtest.in
	@ln -s $(G4INSTALL)/examples/novice/N04/exampleN04.EMtest.out \
          test104/test104.EMtest.out
	@cd tools/lib; $(MAKE)
	@for dir in $(HADEXERDIRS); \
	  do (cd $$dir; $(MAKE) -f GNUmakefile.hadronic.exerciser); done
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done; \
	  for example in $(EXAMPLES); do \
	  testdir=test1$$example; \
	  testin=test1$$example.in; \
	  testout=test1$$example.out; \
	  exdir=$(G4INSTALL)/examples/novice/N$$example; \
	  exin=exampleN$$example.in; \
	  exout=exampleN$$example.out; \
	  rm -f $$testdir/$$testin; \
	  rm -f $$testdir/$$testout; \
	  rm -f $$testdir; \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test1$$example; \
	  ln -s $$exdir $$testdir; \
	  ln -s $$exdir/$$exin $$testdir/$$testin; \
	  ln -s $$exdir/$$exout $$testdir/$$testout; \
	  ln -s $(G4WORKDIR)/bin/$(G4SYSTEM)/exampleN$$example \
	    $(G4WORKDIR)/bin/$(G4SYSTEM)/test1$$example; \
	  (cd $$testdir; $(MAKE)); \
	done

clean:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean); done
	@rm -f test104/test104.EMtest.in
	@rm -f test104/test104.EMtest.out
	@for example in $(EXAMPLES); do \
	  testdir=test1$$example; \
	  exdir=$(G4INSTALL)/examples/novice/N$$example; \
	  rm -f $$testdir/$$testin; \
	  rm -f $$testdir/$$testout; \
	  rm -f $$testdir; \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test1$$example; \
	  ln -s $$exdir $$testdir; \
	  ln -s $$exdir/$$exin $$testdir/$$testin; \
	  ln -s $$exdir/$$exout $$testdir/$$testout; \
	  ln -s $(G4WORKDIR)/bin/$(G4SYSTEM)/exampleN$$example \
	  (cd $$testdir; $(MAKE) clean); \
	  rm -f $$testdir/$$testin; \
	  rm -f $$testdir/$$testout; \
	  rm -f $$testdir; \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test1$$example; \
	done

clean_libs:
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE) clean_libs); done
	@rm -f test104/test104.EMtest.in
	@rm -f test104/test104.EMtest.out
	@for example in $(EXAMPLES); do \
	  testdir=test1$$example; \
	  exdir=$(G4INSTALL)/examples/novice/N$$example; \
	  rm -f $$testdir/$$testin; \
	  rm -f $$testdir/$$testout; \
	  rm -f $$testdir; \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test1$$example; \
	  ln -s $$exdir $$testdir; \
	  ln -s $$exdir/$$exin $$testdir/$$testin; \
	  ln -s $$exdir/$$exout $$testdir/$$testout; \
	  ln -s $(G4WORKDIR)/bin/$(G4SYSTEM)/exampleN$$example \
	  (cd $$testdir; $(MAKE) clean_libs); \
	  rm -f $$testdir/$$testin; \
	  rm -f $$testdir/$$testout; \
	  rm -f $$testdir; \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test1$$example; \
	done

