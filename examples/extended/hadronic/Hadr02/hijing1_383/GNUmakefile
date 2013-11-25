
#CERNROOT:= /home/khelwagd

ifeq "$(TYPE)" ""
TYPE := $(shell uname)
FC              =	gfortran
LD		=       gfortran
RM              =       rm
O               =       .o
FFLAGS		=       -fPIC -fno-automatic
LDFLAGS		=	"-L/usr/lib -lgfortran"
SYSTEMFILES     =       $(SRCGNU)      
endif

#EXTRALIBS =-L${CERNROOT}/cernlib -lmathlib -lkernlib -lpacklib


# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	Main_HIJING
endif

SRC 		= 	hijing1.383.f hipyset1.35.f\
                        G4HIJING_Block_DATA.f\
			$(MAIN).f

#                	G4HIJING_RANDOM.f\
	   		

# -------------------------------------------------

OBJDIR		=	obj_$(TYPE)
SRCFILES 	= 	$(SRC) $(INC) $(SYSTEMFILES) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC) $(SYSTEMFILES))))
TARGET		=	$(MAIN).$(TYPE)
INSTPATH	=	./

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)
all.F: $(SRC)
	cat $^ > all.F 

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(LD) $(OBJECTS) $(LIBS) $(LDFLAGS) -o $(TARGET)

clean:		
		-rm $(OBJECTS) 

distclean:	
		-rm $(TARGET) 
		-rm -r obj_*
		-rm *~

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)/$(basename $(TARGET))
		strip $(INSTPATH)/$(basename $(TARGET))

# --------------- Dependencies -------------------

