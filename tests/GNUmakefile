# $Id: GNUmakefile,v 1.1 1999-01-08 16:34:11 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

ifndef G4INSTALL
  G4INSTALL = ..
endif
 
include $(G4INSTALL)/config/architecture.gmk

#Genuine tests...
SUBDIRS = test06 test09 test12 test13 test01 test03 test05 test07 test08 test10 test02 test11
 
# test07 test101 test102 test103 test104 test105 test106

# Tests with hadronic exercisers.
HADEXERDIRS = test12 test13 test02 test11 

# Symbolic links to examples...
EXAMPLES = 01 02 03 05 06 04

.PHONY : all clean clean_libs

all:
	@for dir in $(HADEXERDIRS); \
	  do (cd $$dir; $(MAKE) -f GNUmakefile.hadronic.exerciser); done
	@for dir in $(SUBDIRS); do (cd $$dir; $(MAKE)); done
	@for example in $(EXAMPLES); do \
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
