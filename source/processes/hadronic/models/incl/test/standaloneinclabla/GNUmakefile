# $Id: GNUmakefile,v 1.1 2007-09-11 13:31:18 miheikki Exp $
#G4TARGET := cascade
#G4TARGET := benchmarks
G4TARGET := inclAblaStandAlone
#G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/binmake.gmk

# Use ROOT for data analysis
#CPPFLAGS += $(shell root-config --cflags)
#LDFLAGS  += $(shell root-config --libs)

#CPPFLAGS += -DUSEROOT

b1:  # benchmark 1
	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/benchmarks 1 >> b1.out; date; cat b1.out 

b2:  # test run for cascade  parameters: runId #events particle energy [MeV] A Z
#	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 1 1000 1 1000 27 13 > ../data/b2.out; date; head ../data/b2.out; wc ../data/b2.out 
#	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 1 1000 1 3000 27 13 > b2.out; date; head b2.out; 


incl:
	rm -f incl.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/testIncl 1  > incl.out; date; head incl.out; 

# 15 MeV n + 54Fe
testFe54Neutron15MeV:
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclAblaStandAlone 54 26 2 15 10000 Fe54Neutron15MeV.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/out2root Fe54Neutron15MeV.out ../data/testrunFe54Neutron15MeV_withabla.root

# 44 MeV p + 93Nb
testNb93Proton44MeV:
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclAblaStandAlone 93 41 1 44 10000 Nb93Proton44MeV.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/out2root Nb93Proton44MeV.out ../data/testrunNb93Proton44MeV_withabla.root

# 62 p +197Au
testAu197Proton62MeV:
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclAblaStandAlone 197 79 1 62 10000 Au197Proton62MeV.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/out2root Au197Proton62MeV.out ../data/testrunAu197Proton62MeV_withabla.root

# 90 p + 209Bi  (ks. kuva ja sen skripti)
testBi209Proton90MeV:
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclAblaStandAlone 209 83 1 90 10000 Bi209Proton90MeV.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/out2root Bi209Proton90MeV.out ../data/testrunBi209Proton90MeV_withabla.root

# 113 p + 208Pb
testPb208Proton113MeV:
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclAblaStandAlone 208 82 1 113 10000 Pb208Proton113MeV.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/out2root Pb208Proton113MeV.out ../data/testrunPb208Proton113MeV_withabla.root

# 160 p + 90Zr
testZr90Proton160MeV:
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclAblaStandAlone 90 40 1 160 10000 Zr90Proton160MeV.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/out2root Zr90Proton160MeV.out ../data/testrunZr90Proton160MeV_withabla.root

# 260 p +90Zr
testZr90Proton260MeV:
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclAblaStandAlone 90 40 1 260 10000 Zr90Proton260MeV.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/out2root Zr90Proton260MeV.out ../data/testrunZr90Proton260MeV_withabla.root

# 160 p + 208 Pb 
testPb208Proton160MeV:
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclAblaStandAlone 208 82 1 160 10000 Pb208Proton160MeV.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/out2root Pb208Proton160MeV.out ../data/testrunPb208Proton160MeV_withabla.root

# 260  p + 208 Pb
testPb208Proton260MeV:
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclAblaStandAlone 208 82 1 260 10000 Pb208Proton260MeV.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/out2root Pb208Proton260MeV.out ../data/testrunPb208Proton260MeV_withabla.root

# 300 MeV p + 56Fe
testFe56Proton300MeV:
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclAblaStandAlone 56 26 1 300 10000 Fe56Proton300MeV.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/out2root Fe56Proton300MeV.out ../data/testrunFe56Proton300MeV_withabla.root

# 1 GeV p + 208Pb
testPb208Proton1GeV:
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclAblaStandAlone 208 82 1 1000 10000 Pb208Proton1GeV.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/out2root Pb208Proton1GeV.out ../data/testrunPb208Proton1GeV_withabla.root

# 425 n + 63Cu
testCu63Neutron425MeV:
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclAblaStandAlone 63 29 2 425 10000 Cu63Neutron425MeV.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/out2root Cu63Neutron425MeV.out ../data/testrunCu63Neutron425MeV_withabla.root

# 1 GeV p + 208Pb
testevap:
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclAblaStandAlone 208 82 1 1000 10000 testevap.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/out2root testevap.out testevap.root
	root -l compare.C
	gv testevap.ps &

testrun:
#	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclStandAlone 208 82 1 1000 1 test.out > console.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/inclStandAlone 27 13 1 1000 1 test.out > console.out

alltests: testFe54Neutron15MeV testNb93Proton44MeV testAu197Proton62MeV testBi209Proton90MeV \
	testPb208Proton113MeV testZr90Proton160MeV testZr90Proton260MeV testPb208Proton160MeV \
	testPb208Proton260MeV testFe56Proton300MeV testCu63Neutron425MeV testPb208Proton1GeV

help:
	@echo "First build the test driver by issuing command: make clean && make"
	@echo "Make sure that the program out2root (in directory test/out2root) is compiled."
	@echo "--------------------------------------------------------------------------------"
	@echo "Available test runs:"
	@echo "alltests: run all tests below"
	@echo "testFe54Neutron15MeV: 15 MeV n + 54Fe, 10000 events"
	@echo "testNb93Proton44MeV: 44 MeV p + 90Nb, 10000 events"
	@echo "testAu197Proton62MeV: 62 MeV p + 197Au, 10000 events"
	@echo "testBi209Proton90MeV: 90 MeV p + 209Bi, 10000 events"
	@echo "testPb208Proton113MeV: 113 MeV p + 208Pb, 10000 events"
	@echo "testZr90Proton160MeV: 160 MeV p + 90Zr, 10000 events"
	@echo "testZr90Proton260MeV: 260 MeV p + 90Zr, 10000 events"
	@echo "testPb208Proton160MeV: 160 MeV p + 208Pb, 10000 events"
	@echo "testPb208Proton260MeV: 260 MeV p + 208Pb, 10000 events"
	@echo "testFe56Proton300MeV: 300 MeV p + 56Fe, 10000 events"
	@echo "testCu63Neutron425MeV: 425 MeV n + 63Cu, 10000 events"
	@echo "testPb208Proton1GeV: 1GeV p + 208Pb, 10000 events"


