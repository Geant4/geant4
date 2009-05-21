name := mainStatAccepTest
G4TARGET := $(name)
G4EXLIB := true

include $(G4INSTALL)/config/architecture.gmk

.PHONY: all
 
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

ifdef G4ANALYSIS_USE  
  CPPFLAGS += `aida-config --include`
  LDFLAGS  += `aida-config --lib`
  LOADLIBS += `aida-config --lib`
endif

ifdef FULLY_STATIC_EXECUTABLE
  LDFLAGS  += -static
endif 
