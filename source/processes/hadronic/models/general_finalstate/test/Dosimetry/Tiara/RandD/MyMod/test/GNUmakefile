name := NeutronTest1
G4TARGET := $(name)
G4EXLIB := true

.PHONY: all
all: lib bin
include $(G4INSTALL)/config/architecture.gmk

# ifdef __USE_BIASING__
CPPFLAGS += -D__USE_BIASING__
# endif

CPPFLAGS += $(C_FLAGS) -Df2cFortran -D__DEBUG__\
		-I$(LHCXX_REL_DIR)/include\
		-I$(CERN_ROOT)/include/cfortran\
		-I$(G4BASE)/processes/hadronic/models/leading_particle/include\
		-I../include

LDFLAGS += -lG4hadronic_leading_particle -lG4hadronic_inc

ifdef G4DEBUG
LDFLAGS += -L/afs/cern.ch/sw/contrib/Mesa/3.2/Linux/lib\
	   -L$(LHCXX_REL_DIR)/lib\
	   -L$(CERN_ROOT)/lib\
	   -L/usr/X11R6/lib\
	   -l$(CLHEP_LIB)\
	   -lAIDA_Annotation\
	   -lAIDA_HBook\
	   -lCHBook\
	   -lAIDA_Plotter\
	   -lVector\
	   -lHepUtilities\
	   -lNeutronTest1\
	   -lqp -lqpview\
	   -lqt-mt\
	   -lgraflib\
	   -lgrafX11\
	   -lpacklib\
	   -lX11\
	   -lnsl\
	   -lcrypt\
	   -lg2c\
	   -ldl
else
LDFLAGS += -L/afs/cern.ch/sw/contrib/Mesa/3.2/Linux/lib\
	   -L$(LHCXX_REL_DIR)/lib\
	   -L$(CERN_ROOT)/lib\
	   -L/usr/X11R6/lib\
	   -l$(CLHEP_LIB)\
	   -lAIDA_Annotation\
	   -lAIDA_HBook\
	   -lCHBook\
	   -lAIDA_Plotter\
	   -lVector\
	   -lHepUtilities\
	   -lNeutronTest1\
	   -lqp -lqpview\
	   -lqt-mt\
	   -lgraflib\
	   -lgrafX11\
	   -lpacklib\
	   -lX11\
	   -lnsl\
	   -lcrypt\
	   -lg2c\
	   -ldl
endif

include $(G4INSTALL)/config/binmake.gmk
