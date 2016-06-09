# $Id: GNUmakefile,v 1.1.1.1 2003/12/01 10:50:22 hpw Exp $
# -----------------------------------------------------------------

SUBDIR1 = FTFC FTFP LBE LHEP LHEP_GN LHEP_HP LHEP_LEAD LHEP_LEAD_HP LHEP_PRECO QGSP_BERT
SUBDIR2 = LHEP_PRECO_HP Packaging QGSC QGSC_LEAD QGSC_LEAD_HP QGSP QGSP_GN QGSP_HP LHEP_BERT LHEP_BIC QGSP_BIC

.PHONY: all clean

all:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE)); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE)); done;:

clean:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE) clean); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) clean); done;:
