# -----------------------------------------------------------------

SUBDIR =  gflash1 gflash2 gflash3 gflasha

.PHONY: all bin clean clean_bin debug

all:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) all); done;:

bin:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) bin); done;:

clean:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) clean); done;:

clean_bin:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) clean_bin); done;:

debug:
	@for dir in $(SUBDIR); do (cd $$dir && $(MAKE) debug); done;:

