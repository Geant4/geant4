# ===========================================================================
#  Makefile urQMD                                    Henning Weber 10.09.1996
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile		[uqmd.f]
##			TYPE	=	operating system	['uname']
##                      LINK    =       linking additional libraries
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install		make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  


# ----------------------------------------- 

ifeq "$(TYPE)" ""
   TYPE := $(shell uname)
endif


ifeq "$(TYPE)" "AIX" 

FC      	= 	xlf
LD              =       xlf
RM              =       rm
O		=	.o
FFLAGS          =       -O5 -qextname -qstrict -qipa=partition=large
LDFLAGS         =       -O5 -qextname -qstrict -qipa=partition=large
SYSTEMFILES	=	$(SRCAIX)

endif


ifeq "$(TYPE)" "Linux" 

FC              =       g77
LD              =       g77
RM              =       rm 
O               =       .o
FFLAGS		=	-g -Wall -Wsurprising -fbounds-check 
LDFLAGS         =       -g
SYSTEMFILES     =       $(SRCGNU)

endif

ifeq "$(TYPE)" "G4INTERFACE"
FC		= 	gfortran
LD		= 	gfortran
RM		=	rm
O		=	.o
FFLAGS		=	-fPIC -c -C -fno-automatic -fdefault-real-8 -ff2c \
			-fno-second-underscore -O -ffixed-line-length-132 \
			-frepack-arrays
LDFLAGS		=	-O
SYSTEMFILES	=	$(SRCGNU)
ifdef G4DEBUG
	FFLAGS  += -g
endif
endif

ifeq "$(TYPE)" "Intel" 

FC              =       ifc -O3 -axMK -tpp6 -w95 -ipo
LD              =       ifc -Vaxlib -O3 -axMK -tpp6 -ipo
RM              =       rm 
O               =       .o
FFLAGS		=	
LDFLAGS         =       
SYSTEMFILES     =       $(SRCINT) $(SRCERF)

endif


ifeq "$(TYPE)" "HP-UX"

FC              =       f77
LD		=	f77
RM		=	rm
O               =       .o
FFLAGS          =       -w
LDFLAGS         =       -g
SYSTEMFILES     =       $(SRCGNU)

endif

ifeq "$(TYPE)" "IRIX64"

FC              =       f77
LD		=	f77
RM		=	rm
O               =       .o
FFLAGS		=	-mips4 -64 -r10000  
LDFLAGS         =       -mips4 -64 -r10000
SYSTEMFILES     =       $(SRCSGI)

endif

ifeq "$(TYPE)" "OSF1"

FC              =       f77
LD              =       f77
RM              =       rm 
O               =       .o 
FFLAGS          =        -C -align dcommons  
LDFLAGS         =        -C -align dcommons 
SYSTEMFILES     =       $(SRCALPHA)

endif

ifeq "$(TYPE)" "SunOS"
 
FC          	=       f77
LD 		=       f77
RM  		=       rm
O  	 	=       .o
FFLAGS      	=       -w
LDFLAGS    	=       -g
SYSTEMFILES  	=       $(SRCSUN)
 
endif


# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	UrqmdMain
endif

SRC 		= 	addpart.f angdis.f anndec.f blockres.f cascinit.f \
	   		coload.f dectim.f delpart.f detbal.f dwidth.f \
	   		getmass.f getspin.f init.f input.f ityp2pdg.f error.f\
	   		make22.f numrec.f output.f paulibl.f jdecay2.f iso.f\
	   		proppot.f saveinfo.f scatter.f siglookup.f string.f \
	   		tabinit.f whichres.f boxprg.f urqmd.f UrqmdLink.f $(MAIN).f\
                        g4urqmdblockdata.f
INC		= 	colltab.f comnorm.f coms.f comres.f comstr.f comwid.f \
			freezeout.f inputs.f newpart.f options.f boxinc.f \
			outcom.f 
#SRCAIX		=	ri6000ranf.f risctime.f
SRCAIX		=	genranf.f
SRCGNU		=	gnuranf.f
SRCALPHA	=	alpharanf.f
SRCSGI		=	genranf.f
SRCSUN		=	gnuranf.f
SRCINT		=	intranf.f
SRCERF		=	erf.f

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
		-rm $(OBJECTS) loadmap.out

distclean:	
		-rm $(TARGET) loadmap.out
		-rm -r obj_*
		-rm *~

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)/$(basename $(TARGET))
		strip $(INSTPATH)/$(basename $(TARGET))

# --------------- Dependencies -------------------

./addpart.f:	coms.f newpart.f
./angdis.f:	options.f
./anndec.f: 	comres.f comwid.f options.f
./blockres.f: 	comres.f
./boxprg.f:     boxinc.f coms.f comres.f options.f
./cascinit.f: 	coms.f inputs.f options.f 
./coload.f: 	colltab.f coms.f comres.f newpart.f options.f
./dectim.f: 	coms.f colltab.f options.f 
./delpart.f:	coms.f comres.f newpart.f
./detbal.f: 	coms.f comres.f comwid.f newpart.f options.f
./dwidth.f:	coms.f comres.f comwid.f options.f
./getmass.f:    comnorm.f comres.f
./getspin.f:	comres.f 
./init.f: 	boxinc.f colltab.f coms.f comres.f freezeout.f inputs.f options.f
./input.f: 	boxinc.f coms.f comres.f inputs.f options.f
./ityp2pdg.f:	comres.f
./iso.f:	comres.f newpart.f
./jdecay2.f:	newpart.f
./make22.f:	comres.f comwid.f newpart.f options.f
./output.f: 	boxinc.f coms.f comres.f freezeout.f inputs.f newpart.f options.f
./paulibl.f: 	coms.f options.f
./proppot.f: 	boxinc.f coms.f options.f
./saveinfo.f: 	coms.f
./scatter.f: 	coms.f comres.f freezeout.f newpart.f options.f
./siglookup.f:	comres.f
./string.f:	comres.f comstr.f options.f
./tabinit.f:	comres.f comwid.f
./urqmd.f: 	boxinc.f colltab.f coms.f comres.f inputs.f newpart.f options.f
./whichres.f: 	comres.f options.f
