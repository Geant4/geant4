all:	# :::
	grep ":::" GNUmakefile

man:	# ::: Make manual.
	make clean; latex manCascade; bibtex manCascade; latex manCascade; latex manCascade; dvips manCascade.dvi -o manCascade.ps; gv manCascade.ps

pub:	# ::: Make publication.
	make clean; latex pubCascade; bibtex pubCascade; latex pubCascade; latex pubCascade; dvips pubCascade.dvi -o pubCascade.ps; gv pubCascade.ps


clean:	# ::: clean directory
	rm -f \#*.*\# \#*\# \#* .*~ .#* *.*~ *~ *.o; rm -f core tmp.txt *.blg *.bbl *.toc *.log *.aux *.dvi