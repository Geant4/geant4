name := hadronics
G4TARGET := $(name)
G4EXLIB := true


ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk


ifdef G4ANALYSIS_USE  
 CPPFLAGS += `aida-config --include`
 LDFLAGS  += `aida-config --lib`
 LOADLIBS += `aida-config --lib` 
endif

export G4ANALYSIS_USE_ROOT  

ifdef G4ANALYSIS_USE_ROOT  
 CPPFLAGS += -fPIC -pthread -I$(ROOTSYS)/include
 LDLIBS  += -L$(ROOTSYS)/lib -lCore -lCint -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -pthread -lm -ldl -rdynamic
endif

