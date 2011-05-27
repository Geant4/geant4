# $Id: GNUmakefile,v 1.4 2008/06/25 17:20:04 kaitanie Exp $
#G4TARGET := cascade
#G4TARGET := benchmarks
G4TARGET := ex1
#G4EXLIB := true
TESTPROGRAM = 
ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/binmake.gmk

TESTPROGRAM = $(G4WORKDIR)/bin/$(G4SYSTEM)/inclAblaStandAlone

# Use ROOT for data analysis
CPPFLAGS += $(shell root-config --cflags)
LDFLAGS  += $(shell root-config --libs)

CPPFLAGS += -DUSEROOT

t0:   #// 197Au79 tested for 30-3000 OK (likewise 54Be26 9Be4 ans 238U92)
	rm -f c.out; date;
	$(G4WORKDIR)/bin/$(G4SYSTEM)/ex1 1 100 1 1200 208 82;

t1:   #// 197Au79 tested for 30-3000 OK (likewise 54Be26 9Be4 ans 238U92)
	mkdir -p data
	rm -f c.out; date;
	$(G4WORKDIR)/bin/$(G4SYSTEM)/ex1 1 10000 1 1200 208 82;

t2:   #// 197Au79 tested for 30-3000 OK (likewise 54Be26 9Be4 ans 238U92)
	mkdir -p data
	rm -f c.out; date;
	$(G4WORKDIR)/bin/$(G4SYSTEM)/ex1 1 10000 1 1200 208 82; # // runId 1 uses INCL/ABLA
	mv data/data.out data/datatmp1.out
	$(G4WORKDIR)/bin/$(G4SYSTEM)/ex1 2 10000 1 1200 208 82; # // runId 2 uses Bertini
	mv data/data.out data/datatmp2.out
	cat data/datatmp1.out > data/data.out
	cat data/datatmp2.out | sed "1d;$d" >> data/data.out  # // Add Bertini but remove header. 
	head data/data.out; wc data/data.out; ls -lh data; ls -lat data; date

vg:
	mkdir data
	valgrind $(G4WORKDIR)/bin/$(G4SYSTEM)/ex1 1 10 1 1200 208 82; 2> vgtest.log

c:
	rm *~
	rm *.log

c9:
	root coulomb.C 

plot1:
	asy plot1.asy; 
