# $Id: GNUmakefile,v 1.3 1999-06-02 08:10:28 stesting Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

ifndef G4INSTALL
  G4INSTALL = ..
endif
 
include $(G4INSTALL)/config/architecture.gmk

#Genuine tests...
SUBDIRS = test01 ## test02 test03
# SUBDIRS += test04 ## Persistency.
SUBDIRS += test05 test06 test07 test08 test09 test10
SUBDIRS += test11 test12
# SUBDIRS += test13 ## Kinetic model
SUBDIRS += test14 ## LowEnegery em.

# Tests with hadronic exercisers.
HADEXERDIRS = test02 test11 test12 test13

# Symbolic links to examples...
EXAMPLES = 01 02 03 04 05 06

.PHONY : all clean clean_libs

all:
	@rm -f test104/test104.EMtest.in
	@rm -f test104/test104.EMtest.out
	@ln -s $(G4INSTALL)/examples/novice/N04/exampleN04.EMtest.in \
          test104/test104.EMtest.in
	@ln -s $(G4INSTALL)/examples/novice/N04/exampleN04.EMtest.out \
          test104/test104.EMtest.out
	@for dir in $(HADEXERDIRS); \
	  do (cd $$dir; $(MAKE) -f GNUmakefile.hadronic.exerciser); done; \
	  for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done; \
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
	@rm test104/test104.EMtest.in
	@rm test104/test104.EMtest.out
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
	@rm test104/test104.EMtest.in
	@rm test104/test104.EMtest.out
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

