# $Id: GNUmakefile,v 1.46 2007-05-20 16:58:01 miheikki Exp $
# $Log: not supported by cvs2svn $
# Revision 1.45  2007/05/20 16:54:36  miheikki
# minor cleaning
#
# Revision 1.44  2007/05/20 15:23:11  miheikki
# added test for G4InuclEvaporation
#
# Revision 1.43  2007/05/20 11:21:56  miheikki
# adding keywords
#
G4TARGET := $(name) # possible targets: cascade benchmarks
#G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/binmake.gmk

b1:  # benchmark 1
	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/benchmarks 1 >> b1.out; date; cat b1.out 

b2:  # test run for cascade  parameters: runId #events particle energy [MeV] A Z
#	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 1 1000 1 1000 27 13 > ../data/b2.out; date; head ../data/b2.out; wc ../data/b2.out 
	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 1 1000 1 3000 27 13 > b2.out; date; head b2.out; 



