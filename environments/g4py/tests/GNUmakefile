# $Id: GNUmakefile,v 1.2 2006-04-25 07:54:28 kmura Exp $
# $Name: not supported by cvs2svn $
# ===========================================================
#   Makefile for building Geant4Py modules
# ===========================================================
SUBDIR1 := test00 test07 test08 test09 test10 test11 
SUBDIR1 += test12 test13
SUBDIR2 := test01 test02 test04 test05 test06

.PHONY: all clean

all:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE)); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) all); done;:

clean:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE) clean); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) cleanall); done;:

