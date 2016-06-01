#---------------------------------------------------------------------
#
# file: makefile
# author: XVT-Architect, XVT Software
#
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# System settings

# Modify the following according to your system setup:
#
# XVT_DSP_DIR	path to the xvt directrorie
# CC            the name of your C++ compiler
# LD            linker
# UIL           uil compiler
# SFX           suffix of C++ source files

CC   		= CC
LD   		= CC
UIL  		= /usr/bin/X11/uil
SFX 		= .cxx

#---------------------------------------------------------------------
# Includes:

INCLUDE =  -I. -I./factory \
			-I$(XVT_DSP_DIR)/include \
			-I/usr/include/Motif1.2 \
			-I/usr/include/X11R5

#---------------------------------------------------------------------
# Libraries:

PWRLIBS     = $(XVT_DSP_DIR)/lib/libxvtxmpwr450.a
PTKLIBS     = $(XVT_DSP_DIR)/lib/libxvtxmapi457.a -lxvtxmhb457 -lxvtxmba457
LIBS        = $(PWRLIBS) $(PTKLIBS) -lxvtxmrw -lMrm -lXm -lXt -lX11 -lPW -lm

#---------------------------------------------------------------------
# Flags:

PWRFLAGS    = -DSTDH -D_HPUX_SOURCE -DWSMTF -D__hpux -D__unix -D__STDC__
CCFLAGS 	= -Aa +a1 $(PWRFLAGS) $(INCLUDE)
LDFLAGS     = -Wl,+s +a1 -L$(XVT_DSP_DIR)/lib -L/usr/lib/Motif1.2 -L/usr/lib/X11R5
CURLFLAGS   = -r mtf12 -s $(INCLUDE)
UILFLAGS    = -I/usr/include/Motif1.2

#-------------------------------------------------------------------
# Default make instructions:
     
.SUFFIXES:	.url .uil .uid .o $(SFX) .h 

.url.uil: 	
	$(XVT_DSP_DIR)/bin/curl $(CURLFLAGS) $<

.uil.uid: 	
	$(UIL) $(UILFLAGS) -o $@ $<

$(SFX).o:	
	$(CC) $(CCFLAGS) -o $@ -c $<

#-------------------------------------------------------------------
# Targets for shell:
APPNAME = GUI

# Shell and factory objects
GENERATEDOBJECTS = \
g4xvtgui.o \
g4xvtdoc.o \
XResMgr.o \
g4xvtwin.o \
cstartup.o \
factory/factory.o \
factory/faccb0.o \
factory/faccb1.o \
factory/faccb2.o \
factory/faccb3.o \
factory/faccb4.o \
factory/classes.o 
 
# other application objects
#   List all other object files which are needed for your application here
APPNAMEOBJECTS = 

OBJECTS = $(GENERATEDOBJECTS) $(APPNAMEOBJECTS)
 
$(APPNAME):  $(OBJECTS) $(APPNAME).uil $(APPNAME).uid
	$(LD) $(LDFLAGS) -o $(APPNAME) $(OBJECTS) $(LIBS)

$(APPNAME).uil: $(APPNAME).url factory/images.url factory/strings.url factory/factory.url
