include $(G4INSTALL)/config/sys/$(G4SYSTEM).gmk

%.o: %.f
	@echo "Compiling $< with $(FC) ..."
	$(FC) -c $(FCFLAGS) -o $@ $< 

initpydata.o : initpydata.f

.PHONY: clean

clean:
	@rm initpydata.o

