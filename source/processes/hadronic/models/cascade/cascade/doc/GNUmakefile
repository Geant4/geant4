# Preparing HETC cascade documentation
all:
	# make web pages
	texi2html cascade
	# make pdf document
	#	texi2pdf sample #needs pdftex.fmt file
	# make ps from texinfo file
	tex cascade
	dvips cascade.dvi
	gv cascade.ps &
clean:
	rm *.html *.ps *.dvi *.cp *.cps *.fn *.info *.ky *.log *.pg *.toc *.tp *.vr *.aux	