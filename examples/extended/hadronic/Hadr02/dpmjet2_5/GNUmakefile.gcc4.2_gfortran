#
# makefile for DPMJET25
#
#-frecord-marker=4
#
# LINUX (Intel PC)
  OPT1 = -c -C -fno-automatic -fdefault-real-8 -ff2c \
         -fno-second-underscore -O -ffixed-line-length-132 -frepack-arrays 
  OPT2 = -c -C -fno-automatic -Wuninitialized -fdefault-real-8 -ff2c \
         -fno-second-underscore -O -ffixed-line-length-132 -frepack-arrays 

ifdef G4DEBUG
  OPT1 += -g 
  OPT2 += -g 
endif

# LIB =
#
# fortran compiler, linker/loader:
#
FC = gfortran
#
  DPMJET2_5: dpm25nuc1g4.o dpm25nuc2.o   dpm25nuc3dpm.o dpm25nuc4.o \
             dpm25nuc5.o   dpm25nuc6g4.o dpm25nuc7g4.o  grv98lo.o  \
             dpm25hadri.o  dpm25nulib.o  dpm25histg4.o \
             dpm25diff.o   dpm25eva.o    pythia61.o     dpm25pomg4.o \
             dpm25lap.o    dpm25nonu.o   g4dpmjet_initialise_block_data.o \
	     g4dpmjet_open_nuclear_bin.o g4dpmjet_close_nuclear_bin.o \
	     g4dpmjet_open_fort6.o       g4dpmjet_close_fort6.o 

dpm25nuc1g4.o : dpm25nuc1g4.f
	$(FC) $(OPT1) dpm25nuc1g4.f

dpm25nuc2.o : dpm25nuc2.f
	$(FC) $(OPT1) dpm25nuc2.f

dpm25nuc3dpm.o : dpm25nuc3dpm.f
	$(FC) $(OPT1) dpm25nuc3dpm.f

dpm25nuc4.o : dpm25nuc4.f
	$(FC) $(OPT1) dpm25nuc4.f

dpm25nuc5.o : dpm25nuc5.f
	$(FC) $(OPT1) dpm25nuc5.f

dpm25nuc6g4.o : dpm25nuc6g4.f
	$(FC) $(OPT1) dpm25nuc6g4.f

dpm25nuc7g4.o : dpm25nuc7g4.f
	$(FC) $(OPT1) dpm25nuc7g4.f

grv98lo.o : grv98lo.f
	$(FC) $(OPT1) grv98lo.f

dpm25hadri.o : dpm25hadri.f
	$(FC) $(OPT1) dpm25hadri.f

dpm25nulib.o : dpm25nulib.f
	$(FC) $(OPT1) dpm25nulib.f

dpm25histg4.o : dpm25histg4.f
	$(FC) $(OPT1) dpm25histg4.f

dpm25diff.o : dpm25diff.f
	$(FC) $(OPT1) dpm25diff.f

dpm25eva.o : dpm25eva.f
	$(FC) $(OPT1) dpm25eva.f

pythia61.o : pythia61.f
	$(FC) $(OPT1) pythia61.f

dpm25pomg4.o : dpm25pomg4.f
	$(FC) $(OPT1) dpm25pomg4.f

dpm25lap.o : dpm25lap.f
	$(FC) $(OPT1) dpm25lap.f

dpm25nonu.o : dpm25nonu.f
	$(FC) $(OPT1) dpm25nonu.f

g4dpmjet_initialise_block_data.o : g4dpmjet_initialise_block_data.f
	$(FC) $(OPT2) g4dpmjet_initialise_block_data.f

g4dpmjet_open_nuclear_bin.o : g4dpmjet_open_nuclear_bin.f
	$(FC) $(OPT2) g4dpmjet_open_nuclear_bin.f

g4dpmjet_close_nuclear_bin.o : g4dpmjet_close_nuclear_bin.f
	$(FC) $(OPT2) g4dpmjet_close_nuclear_bin.f

g4dpmjet_open_fort6.o : g4dpmjet_open_fort6.f
	$(FC) $(OPT2) g4dpmjet_open_fort6.f

g4dpmjet_close_fort6.o : g4dpmjet_close_fort6.f
	$(FC) $(OPT2) g4dpmjet_close_fort6.f


