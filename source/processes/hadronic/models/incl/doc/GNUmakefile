# GNUmakefile for building the reference documentation with Doxygen.

VERSION := $(shell git describe --tags)

all: physref

rl:
	echo "\begin{verbatim}" > releaselog.tex
	echo "- INCL/ABLA changes since incl4.2-abla3-rel1:" >> releaselog.tex
	git --no-pager log --no-merges --pretty="format:  o %s" \
	incl4.2-abla3-rel1..HEAD -- ../include/ ../src/ >> releaselog.tex
	echo "" >> releaselog.tex
	echo "\end{verbatim}" >> releaselog.tex

clean:
	rm -f \#*.*\# \#*\# \#* .*~ .#* *.*~ *~ *.o; rm -f core tmp.txt *.blg *.bbl *.toc *.log *.aux *.dvi
	mkdir -p reference/*
	rm -rf reference/*

docball:
	@echo "Building documentation tarball..."
	tar --no-recursion -czf incl-phys-ref.tar.gz GNUmakefile useIncl.tex Incl.tex *.eps

physref: useIncl.tex Incl.tex
	make clean; latex useIncl; bibtex useIncl; latex useIncl; latex useIncl; dvips useIncl.dvi -o useIncl.ps;

referencedocs: clean
	mkdir -p reference/
	cat Doxyfile.in | sed s/INCLABLAVERSION/$(VERSION)/g > Doxyfile
	doxygen

uploaddocs:
	chmod -R a+r reference/*
	chmod a+x reference/*
	chmod a+x reference/html/*
	scp -r reference/html/* kaitanie@shell.dx.fi:~/public_html/incl/
