all:	# :::
	grep ":::" GNUmakefile

man:	# ::: make manual
	latex manCascade.tex; latex manCascade.tex; dvips manCascade.dvi -o manCascade.ps; gv manCascade.ps

pub:	# ::: make publication/paper
	latex pubCascade.tex; latex pubCascade.tex; dvips pubCascade.dvi -o pubCascade.ps; gv pubCascade.ps
