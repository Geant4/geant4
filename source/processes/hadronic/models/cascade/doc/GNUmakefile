all:	# :::
	grep ":::" GNUmakefile

man:	# ::: Make manual.
	make clean; latex manCascade; bibtex manCascade; latex manCascade; latex manCascade; dvips manCascade.dvi -o manCascade.ps;

pub:	# ::: Make publication.
	make clean; latex pubCascade; bibtex pubCascade; latex pubCascade; latex pubCascade; dvips pubCascade.dvi -o pubCascade.ps; 


clean:	# ::: clean directory
	rm -f \#*.*\# \#*\# \#* .*~ .#* *.*~ *~ *.o; rm -f core tmp.txt *.blg *.bbl *.toc *.log *.aux *.dvi
physRef:
	make clean; latex useCascade; bibtex useCascade; latex useCascade; latex useCascade; dvips useCascade.dvi -o useCascade.ps; 