# $Id: GNUmakefile,v 1.59 2004-01-29 16:10:12 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

ifndef G4INSTALL
  G4INSTALL = ..
endif
 
include $(G4INSTALL)/config/architecture.gmk

# Symbolic links to novice examples using test101, etc...
EXAMPLES = 01 02 03 04 05 06 07

# Symbolic links to extended persistency examples using test401, etc...
#temp PERSISTENCY += 01 02

# Symbolic links to extended electromagnetic examples using test501, etc...
ELECTROMAGNETIC += 1 2 3 4 5 ## Awaiting description.
ELECTROMAGNETIC += 6 7 ## G4GammaConversionToMuons process.
ELECTROMAGNETIC += 8 9 10 ## Transistion radiation and PAI (photoabsorption
                        ## ionisation) model.

# Symbolic links to extended g3tog4 examples using test601, etc...
#temp G3TOG4 += 01 02   ## 01 - clGeometry  02 cltog4 (actually hand coded for now)

# Symbolic links to extended field examples using test701, etc...
FIELD += 01 02 03

# Symbolic links to extended gdml examples using test901, etc...
#temp GDML += 01

# Symbolic links to extended geometry OLAP examples using test1001, etc...
#temp OLAP += 01

# Symbolic links to extended BIAS examples using test801, etc...
BIAS += 01 
#temp BIAS_AIDA += 02
#temp BIAS_PYTHON += 03

#02 03 04 05 06 07


# System Tests...
SUBDIRS += test01
SUBDIRS += test02
# SUBDIRS += test03 ## Optical Processes.  (Derived originally from ExN06.)
                    ## Archived (cvs rm) - was indentical to ExN06 and no
                    ## longer maintained.
# SUBDIRS += test04 ## Persistency.  Retired.
# SUBDIRS += test05
SUBDIRS += test06
SUBDIRS += test07
# SUBDIRS += test08 ## Old em test.  Archived (cvs rm) for now.
                    ## (Action? Replace with transition radiation
                    ## and PAI (photoabsorption ionisation) model.)  (Derived
                    ## from TestEm8.)
SUBDIRS += test09
SUBDIRS += test10
SUBDIRS += test11
SUBDIRS += test12
SUBDIRS += test13
SUBDIRS += test14
SUBDIRS += test15
SUBDIRS += test16
SUBDIRS += test17
SUBDIRS += test18
SUBDIRS += test20
SUBDIRS += test21
SUBDIRS += test22
# SUBDIRS += test23 ## Old test on 'cuts per material', now obsolete 
SUBDIRS += test24
SUBDIRS += test25
SUBDIRS += test26
SUBDIRS += test27
# SUBDIRS += test30 ## Not in testing yet!
# SUBDIRS += test31 ## Not in testing yet!
SUBDIRS += test32
SUBDIRS += test33
SUBDIRS += test40

SUBDIRS_AIDA += test50

# Tests with hadronic exercisers.
HADEXERDIRS = test02 test11 test12 test13 test15 test16 test21 

#temp ADV_EXAMPLES += 01 02  03 04 05 06 07

.PHONY : all examples_links clean clean_libs clean_examples_links

all:
# Making links in the source tree can be a problem if more than one
# jobs try to do it at the same time.  So disable and always do by
# hand with: gmake examples_links_source
###################################################################

# Novice
	@$(ECHO) Making examples...
	@for nn in $(EXAMPLES); do \
	  (cd test1$$nn && $(MAKE)); \
	done;:

# EM
	@$(ECHO) Making EM Tests...
	@for nn in $(ELECTROMAGNETIC); do \
	  (cd test50$$nn && $(MAKE)); \
	done;:

# Persistancy
ifdef G4USE_HEPODBMS
	@$(ECHO) Making Persistency Tests...
	@for nn in $(PERSISTENCY); do \
	  (cd test4$$nn && $(MAKE)); \
	done;:
endif

#temp 
# G3TOG4
#	@$(ECHO) Making G3toG4 Tests...
#	@for nn in $(G3TOG4); do \
#	  (cd test6$$nn && $(MAKE)); \
#	done;:

# Field
	@$(ECHO) Making Field Tests...
	@for nn in $(FIELD); do \
	  (cd test7$$nn && $(MAKE)); \
	done;:

#ifneq (,$(findstring 7.3_2.95.2,$(DEBOPT)))
#
#      @$(ECHO) Making Advanced examples...
#      @for nn in $(ADV_EXAMPLES); do \
#        (cd test10$$nn && $(MAKE)); \
#      done;:
#
#      @$(ECHO) Making AIDA examples...
#      @for nn in $(SUBDIRS_AIDA); do \
#        (cd test$$nn && $(MAKE)); \
#      done;:
#
#      @$(ECHO) Making BIAS AIDA examples...
#      @for nn in $(BIAS_AIDA); do \
#        (cd test8$$nn && $(MAKE)); \
#      done;:
#
#endif

#temp
# GDML
#	@$(ECHO) Making GDML Tests...
#	@for nn in $(GDML); do \
#	  (cd test9$$nn && $(MAKE) XERCESCROOT=$(XERCESCROOT)); \
#	done;:

#temp
# OLAP
#	@$(ECHO) Making OLAP Tests...
#	@for nn in $(OLAP); do \
#	  (cd test10$$nn && $(MAKE)); \
#	done;:

# Biasing
	@$(ECHO) Making Biasing Tests...
	@for nn in $(BIAS); do \
	  (cd test8$$nn && $(MAKE)); \
	done;:

# We don't need it anymore (?)
#	@$(ECHO) Making libG4TestTool.a for System Tests...
#	@cd tools/lib && $(MAKE);:

# Hadronic exercises
	@$(ECHO) Making hadronic exercisers for System Tests...
	@for dir in $(HADEXERDIRS); \
	  do (cd $$dir && $(MAKE) -f GNUmakefile.hadronic.exerciser); done;:

# Tests
	@$(ECHO) Making System Tests...
	@for dir in $(SUBDIRS); do (cd $$dir && $(MAKE)); done;:

# BinLinks
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

################################################################################
examples_links_source:

# Clean source links
	@$(MAKE) clean_examples_links_source

# Novice
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

# EM
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

# Persistency
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

# G3TOG4
	@testdir=test601; \
         exdir=$(G4INSTALL)/examples/extended/g3tog4/clGeometry; \
         ln -s $$exdir $$testdir; \
         ln -s clGeometry.in  $$testdir/test601.in;  \
         ln -s clGeometry.out $$testdir/test601.out; \
         true

	@testdir=test602; \
	exdir=$(G4INSTALL)/examples/extended/g3tog4/cltog4; \
         ln -s $$exdir $$testdir; \
         ln -s cltog4.in  $$testdir/test602.in;  \
         ln -s cltog4.out $$testdir/test602.out; :


# GDML
	@testdir=test901; \
	exdir=$(G4INSTALL)/examples/extended/gdml; \
         ln -s $$exdir $$testdir; :

# OLAP
	@testdir=test1001; \
	exdir=$(G4INSTALL)/examples/extended/geometry/olap; \
         ln -s $$exdir $$testdir; :

################################################################
# Some "truc" to use two input files at the same test directory
################################################################
	@testdir=test34; \
	exdir=$(G4INSTALL)/tests/test33; \
         ln -s $$exdir $$testdir; \
         ln -s MassGeo_TimedApp.mac  $$testdir/test33.in;  \
         ln -s ParallelGeo_TimedApp.mac $$testdir/test34.in; \
         ln -s MassGeo_TimedApp.out  $$testdir/test33.out;  \
         ln -s ParallelGeo_TimedApp.out $$testdir/test34.out; :

################################################################

# Field
	@for nn in $(FIELD); do \
	  testdir=test7$$nn; \
	  exdir=$(G4INSTALL)/examples/extended/field/field$$nn; \
	  ln -s $$exdir $$testdir; \
	  for i in $$exdir/*.in; do \
	    exin=`basename $$i`; \
	    exout=`basename $$i .in`.out; \
	    testin=`echo $$exin | sed s/field/test7/`; \
	    testout=`echo $$exout | sed s/field/test7/`; \
	    ln -s $$exin $$testdir/$$testin; \
	    ln -s $$exout $$testdir/$$testout; \
	  done; \
	done;:

# BIAS
	@for nn in $(BIAS); do \
	  testdir=test8$$nn; \
	  exdir=$(G4INSTALL)/examples/extended/biasing/B$$nn; \
	  ln -s $$exdir $$testdir; \
	  for i in $$exdir/*.cc; do \
	    exin=`basename $$i .cc`.in; \
	    exout=`basename $$i .cc`.out; \
	    testin=`echo $$exin | sed s/exampleB/test8/`; \
	    testout=`echo $$exout | sed s/exampleB/test8/`; \
	    ln -s $$exout $$testdir/$$testout; \
	  done; \
	done;:

# BIAS_AIDA
	@for nn in $(BIAS_AIDA); do \
	  testdir=test8$$nn; \
	  exdir=$(G4INSTALL)/examples/extended/biasing/B$$nn; \
	  ln -s $$exdir $$testdir; \
	  for i in $$exdir/*.cc; do \
	    exin=`basename $$i .cc`.in; \
	    exout=`basename $$i .cc`.out; \
	    testin=`echo $$exin | sed s/exampleB/test8/`; \
	    testout=`echo $$exout | sed s/exampleB/test8/`; \
	    ln -s $$exout $$testdir/$$testout; \
	  done; \
	done;:

# BIAS_PYTHON
	@for nn in $(BIAS_PYTHON); do \
	  testdir=test8$$nn; \
	  exdir=$(G4INSTALL)/examples/extended/biasing/B$$nn; \
	  ln -s $$exdir $$testdir; \
	  for i in $$exdir/*.cc; do \
	    exin=`basename $$i .cc`.in; \
	    exout=`basename $$i .cc`.out; \
	    testin=`echo $$exin | sed s/exampleB/test8/`; \
	    testout=`echo $$exout | sed s/exampleB/test8/`; \
	    ln -s $$exout $$testdir/$$testout; \
	  done; \
	done;:

# ADV_EXAMPLES
	@for nn in $(ADV_EXAMPLES); do \
	  testdir=test8$$nn; \
	  exdir=$(G4INSTALL)/examples/extended/biasing/B$$nn; \
	  ln -s $$exdir $$testdir; \
	  for i in $$exdir/*.cc; do \
	    exin=`basename $$i .cc`.in; \
	    exout=`basename $$i .cc`.out; \
	    testin=`echo $$exin | sed s/exampleB/test8/`; \
	    testout=`echo $$exout | sed s/exampleB/test8/`; \
	    ln -s $$exout $$testdir/$$testout; \
	  done; \
	done;:

###########################################################################
examples_links_binaries:

# Clean up!
#	@$(MAKE) clean_examples_links_binaries

# Novice
	@for nn in $(EXAMPLES); do \
	  ln -s exampleN$$nn $(G4WORKDIR)/bin/$(G4SYSTEM)/test1$$nn; \
	done;:

# EM
	@for nn in $(ELECTROMAGNETIC); do \
	  ln -s TestEm$$nn $(G4WORKDIR)/bin/$(G4SYSTEM)/test50$$nn; \
	done;:

#temp
# Persistency
#	@for nn in $(PERSISTENCY); do \
#	  exdir=$(G4INSTALL)/examples/extended/persistency/PersistentEx$$nn; \
#	  ln -s $$exdir/PersistentEx$${nn}run \
#	    $(G4WORKDIR)/bin/$(G4SYSTEM)/test4$$nn; \
#	done;:

# Field
	@for nn in $(FIELD); do \
	  ln -s field$$nn $(G4WORKDIR)/bin/$(G4SYSTEM)/test7$$nn; \
	done;:

# BIAS
	@for nn in $(BIAS); do \
	  ln -s exampleB$$nn $(G4WORKDIR)/bin/$(G4SYSTEM)/test8$$nn; \
	done;:

#temp
# BIAS_AIDA
#	@for nn in $(BIAS_AIDA); do \
#	  ln -s exampleB$$nn $(G4WORKDIR)/bin/$(G4SYSTEM)/test8$$nn; \
#	done;:

#temp
# BIAS_PYTHON
#	@for nn in $(BIAS_PYTHON); do \
#	  ln -s exampleB$$nn $(G4WORKDIR)/bin/$(G4SYSTEM)/test8$$nn; \
#	done;:

#temp
# G3TOG4
#	@ln -s clGeometry $(G4WORKDIR)/bin/$(G4SYSTEM)/test601; \
#	ln -s cltog4    $(G4WORKDIR)/bin/$(G4SYSTEM)/test602; \

#temp
# OLAP
#	ln -s olapex    $(G4WORKDIR)/bin/$(G4SYSTEM)/test1001; \

#temp
# GDML
#	ln -s gogdml    $(G4WORKDIR)/bin/$(G4SYSTEM)/test901; \

# Truc for test33/test34
	ln -s test33    $(G4WORKDIR)/bin/$(G4SYSTEM)/test34; \
        true

###########################################################################
clean_examples_links_source:

# GDML
	@rm -f test901

# OLAP
	@rm -f test1001

# Truc for test33/test34
	@rm -f test34
	@rm -f test33/test33.in
	@rm -f test33/test33.out
	@rm -f test33/test34.in
	@rm -f test33/test34.out

# Unstandard testname
	@rm -f test104/test104.EMtest.in
	@rm -f test104/test104.EMtest.out

# Novice
	@for nn in $(EXAMPLES); do \
	  testdir=test1$$nn; \
	  for i in $$testdir/$$testdir.*; do \
	    rm -f $$i; \
	  done; \
	  rm -f $$testdir; \
	done;:

# EM
	@for nn in $(ELECTROMAGNETIC); do \
	  testdir=test50$$nn; \
	  for i in $$testdir/$$testdir.*; do \
	    rm -f $$i; \
	  done; \
	  rm -f $$testdir; \
	done;:

# FIELD
	@for nn in $(FIELD); do \
	  testdir=test7$$nn; \
	  for i in $$testdir/$$testdir.*; do \
	    rm -f $$i; \
	  done; \
	  rm -f $$testdir; \
	done;:

# BIAS
	@for nn in $(BIAS); do \
	  testdir=test8$$nn; \
	  for i in $$testdir/$$testdir.*; do \
	    rm -f $$i; \
	  done; \
	  rm -f $$testdir; \
	done;:

# BIAS_AIDA
	@for nn in $(BIAS_AIDA); do \
	  testdir=test8$$nn; \
	  for i in $$testdir/$$testdir.*; do \
	    rm -f $$i; \
	  done; \
	  rm -f $$testdir; \
	done;:

# BIAS_PYTHON
	@for nn in $(BIAS_PYTHON); do \
	  testdir=test8$$nn; \
	  for i in $$testdir/$$testdir.*; do \
	    rm -f $$i; \
	  done; \
	  rm -f $$testdir; \
	done;:

# Percictency
	@for nn in $(PERSISTENCY); do \
	  testdir=test4$$nn; \
	  for i in $$testdir/$$testdir.*; do \
	    rm -f $$i; \
	  done; \
	  rm -f $$testdir; \
	done;:

# G3TOG4
	@rm -f test601/test601.in test601/test601.out test601;:
	@rm -f test602/test602.in test602/test602.out test602;:

########################################################################
clean_examples_links_binaries:

# Novice	
	@for nn in $(EXAMPLES); do \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test1$$nn; \
	done;:

# EM
	@for nn in $(ELECTROMAGNETIC); do \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test50$$nn; \
	done;:

# Persistency
	@for nn in $(PERSISTENCY); do \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test4$$nn; \
	done;:

# Field
	@for nn in $(FIELD); do \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test7$$nn; \
	done;:

# BIAS
	@for nn in $(BIAS); do \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test8$$nn; \
	done;:

# BIAS_AIDA
	@for nn in $(BIAS_AIDA); do \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test8$$nn; \
	done;:

# BIAS_PYTHON
	@for nn in $(BIAS_PYTHON); do \
	  rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test8$$nn; \
	done;:

# G3TOG4
	@rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test601 $(G4WORKDIR)/bin/$(G4SYSTEM)/test602;:

# GDML
	@rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test901 ;:

# OLAP
	@rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test1001;:

# Truc for test33/test34
	@rm -f $(G4WORKDIR)/bin/$(G4SYSTEM)/test34;:
