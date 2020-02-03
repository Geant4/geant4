# -----------------------------------------------------------------

SUBDIR =  B4a B4b B4c B4d

.PHONY: bin clean clean_bin debug

bin:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) bin); done;:
 
clean:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) clean); done;:
 
clean_bin:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) clean_bin); done;:
 
debug:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) debug); done;:
 
