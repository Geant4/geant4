name := G4PreCompoundTest
G4TARGET := $(name)


ifndef G4INSTALL
  G4INSTALL = ../../../../../../..
endif

G4EXLIB := true


CPPFLAGS := -DPRECOMPOUND_TEST -I./Utilities/gzio/include
EXTRALIBS := -lz
#EXTRA_LINK_DEPENDENCIES := $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/libgzio.so


.PHONY: all 

all:  lib bin clear 


include $(G4INSTALL)/config/binmake.gmk



#$(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/libgzio.so:
#	cd Utilities/gzio; $(MAKE)



clear: 
	@rm -f ./*~
	@rm -f ./core
	@rm -f include/*~
	@rm -f src/*~

