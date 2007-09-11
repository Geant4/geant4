# GNUmakefile for building the reference documentation with Doxygen.

all: physicsreference

clean:
	mkdir -p reference/*
	rm -rf reference/*

physicsreference: physicsreference.tex
	latex physicsreference
	latex physicsreference
	dvips -o physicsreference.ps physicsreference.dvi
#	pdflatex physicsreference
#	pdflatex physicsreference

referencedocs: clean
	mkdir -p reference/*
	doxygen

uploaddocs:
	chmod -R a+r reference/
	chmod a+x reference/
	chmod a+x reference/html/
	scp -r reference/ kaitanie@melkki.cs.helsinki.fi:~/public_html/inclabla/
