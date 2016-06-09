# $Id: GNUmakefile,v 1.6 2006/06/21 18:09:19 vnivanch Exp $
# -----------------------------------------------------------------

SUBDIR1 = FTFC FTFP LBE LHEP LHEP_LEAD LHEP_LEAD_HP LHEP_PRECO QGSP_BERT QGSP_BERT_HP QGSP_EMV LHEP_EMV QGSP_EMX QBBC
SUBDIR2 = LHEP_PRECO_HP Packaging QGSC QGSC_LEAD QGSC_LEAD_HP QGSP QGSP_HP LHEP_BERT LHEP_BERT_HP LHEP_BIC LHEP_BIC_HP QGSP_BIC

.PHONY: all clean

all:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE)); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE)); done;:

clean:
	@for dir in $(SUBDIR1); do (cd $$dir && $(MAKE) clean); done;:
	@for dir in $(SUBDIR2); do (cd $$dir && $(MAKE) clean); done;:
