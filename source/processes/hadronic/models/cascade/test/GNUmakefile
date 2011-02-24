# $Id: GNUmakefile,v 1.60 2008-09-30 22:05:59 miheikki Exp $
#G4TARGET := cascade
#G4TARGlET := interface
G4TARGET := benchmarks
#G4TARGET := evaporation
#G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/binmake.gmk
# enum particleType { nuclei = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, foton = 10 };

unit:
	(cd unit; gmake;) //AH:::

t1:   #// 197Au79 tested for 30-3000 OK (likewise 54Be26 9Be4 ans 238U92)
	rm -f c.out; date;
	echo "0 20000 1 62 27 13" > c.out  # header info for scripts/analalyzeEvents.C : runID totNumOfCollisions bulletType momZ A Z
	$(G4INSTALL)/bin/$(G4SYSTEM)/interface 1 10000 1 62 27 13 >> c.out; # // runId 1 uses standard interface
	$(G4INSTALL)/bin/$(G4SYSTEM)/interface 2 10000 1 62 27 13 >> c.out; # // runId 0 used development version G4IBertini
#	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade   2 1000 1 1000 197 79  >> c.out; # Fe (worst case 20 MeV p-Fe p)
	mv c.out data/data.out; head data/data.out; wc data/data.out; ls -lh data; ls -lat data; date

i1:
	rm -f l.out; date; 
	echo "1 1000 1 1000 197 79" >l.out  # header info for scripts/analalyzeEvents.C 
	$(G4INSTALL)/bin/$(G4SYSTEM)/interface 1 1000 1 1000 197 79  >> l.out; # [l]og
	mv l.out data/l.out; head data/l.out; wc data/l.out; ls -lh data; ls -lat data; date 


c9:  # coulomb test 
	rm -f c.out; date;
	echo "0 11000 1 40 2 1" >c.out  # header info for scripts/analalyzeEvents.C : runID totNumOfCollisions bulletType momZ A Z
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 1  1000 1 40 2   1  >>  c.out; # D
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 2  1000 1 40 9   4  >> c.out; # Be
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 3  1000 1 40 19  9  >> c.out; # F  (2.9 MeV)
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 4  1000 1 40 27  13 >> c.out; # Al
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 5  1000 1 40 54  26 >> c.out; # Fe (worst case 20 MeV p-Fe p) 6.5MeV Coulob barrier
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 6  1000 1 40 89  39 >> c.out; # Y
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 7  1000 1 40 120 50 >> c.out; # Sn
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 8  1000 1 40 169 69 >> c.out; # Tm 12.7 MeV
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 9  1000 1 40 197 79 >> c.out; # Au (best case 20 MeV p-Au p spectrum ~ ok)
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 10 1000 1 40 209 83 >> c.out; # Bi 
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 11 1000 1 40 238 92 >> c.out; # U (15.3 MeV counted by G4InuclCollider::inelasticInteractionPossible
	mv c.out data/data.out; head data/data.out; wc data/data.out; ls -lh data; ls -lat data; date


c9p:   # c9 [p]lots
	cd scripts; root coulomb.C

e1:  # evaporation test1
	rm -f e1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/testEvaporation  >> e1.out; date

b1:  # benchmark 1
	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/benchmarks 1 >> b1.out; date; cat b1.out 

c8:  # 26-Iron Fe
	rm -f c.out; date; 
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 100 10000 1 100 54 26 > c.out;
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 200 10000 1 200 54 26 >> c.out;
	cp c.out data/data.out; head c.out; date 

# castaneda, precompound ... fig 3. 28 Nickel Ni
c7:  # 26-Iron Fe
	rm -f c.out; date; 
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 29 50000 1 29 54 26 > c.out;
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 62 50000 1 62 54 26 >> c.out;
	cp c.out data/data.out; head c.out; date 

c6:
	rm -f c.out; date; 
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 58 10000 2 60 58 28 > c.out;
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 64 10000 2 64 58 28 >> c.out;
	cp c.out data/data.out; head c.out; date 

c5:  #  M. Mann Phys. Rev C Vol. 28, No.4  1983 pp.1487  fig. 17 50 - Tin - Sn
	rm -f c5.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 25 10000 1 25 120 50 > c5.out;
	$(G4INSTALL)/bin/$(G4SYSTEM)/cascade 35 10000 1 35 120 50 >> c5.out; 
	cp c5.out data/data.out; head c5.out; date 


c4:  # compare with blann p.1657
	rm -f c4.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 90 10000 5 90 197 79 > c4.out; date; head c4.out; cp c4.out data/data.out 

c3:  # compare with kalbach05 Fig.11 c)
	rm -f c3.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 62 10000 1 62 197 79 > c3.out; date; head c3.out; cp c3.out data/data.out 

c2:  # compare with heikkinen07::: 1(b) Bismut Bi isotope 100% neutron production
	rm -f c2.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 90 10000 1 90 209 83 > c2.out; date; head c2.out; cp c2.out data/data.out 

#   Tin Sn
#	rm -f c2.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 35 10000 1 35 119 50 > c2.out; date; head c2.out; cp c2.out data/data.out 

c1:  # test run for cascade  parameters: runId #events particle energy [MeV] A Z
	rm -f b1.out; date; $(G4INSTALL)/bin/$(G4SYSTEM)/cascade 1000 1000 1 1000 27 13 > c1.out; date; head c1.out; cp c1.out data/data.out 

eva:
	make clean; make; make e1


# cd $G4INSTALL/source/materials/include
# g Fe *  >>> G4StaticSandiaData.hh:377:// 26 - Iron - Fe - 7
