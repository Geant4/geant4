all:
	latex manCascade.tex; manCascade.tex; manCascade.tex; dvips manCascade.dvi -o manCascade.ps; gv manCascade.ps

release: all
	ps2pdf manCascade.ps; latex2html manCascade.tex

