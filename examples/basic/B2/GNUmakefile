# $Id: GNUmakefile,v 1.97 2010-09-30 09:05:24 gcosmo Exp $
# -----------------------------------------------------------------

SUBDIR =  B2a B2b

.PHONY: bin clean clean_bin debug

bin:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) bin); done;:
 
clean:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) clean); done;:
 
clean_bin:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) clean_bin); done;:
 
debug:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) debug); done;:
 
