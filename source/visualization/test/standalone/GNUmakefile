# Usage: make name=<program>
# where there exists a file containing main() with name program.cc.
# E.g:   make name=standalone
#        make name=standalone clean

# Check user has supplied name and no target...
ifdef name

G4TARGET := $(name)
G4EXLIB := true

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

else # User has not supplied name...
any $(MAKECMDGOALS):
	@echo "**** Usage: make name=<program>"
	@echo "**** where there exists a file containing main() with name program.cc."
	@echo "**** E.g:   make name=standalone"
	@echo "****        make name=standalone clean"
endif
