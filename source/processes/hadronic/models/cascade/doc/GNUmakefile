all:	# :::
	grep ":::" GNUmakefile

man:	# ::: Make manual.
	latex manCascade.tex; latex manCascade.tex; dvips manCascade.dvi -o manCascade.ps; gv manCascade.ps

pub:	# ::: Make publication.
	latex pubCascade.tex; latex pubCascade.tex; bibtex pubCascade.tex; latex pubCascade.tex; latex pubCascade.tex; dvips pubCascade.dvi -o pubCascade.ps; gv pubCascade.ps
