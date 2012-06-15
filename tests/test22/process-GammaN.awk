#!/usr/bin/awk -f
#
#  Given output of test22 with '/tracking/verbose 1', produce a table of the
#  number of times each interaction occurs.
#  Displays the Material, primary particle and Energy (value, unit), 
#   in addition to output data.
#
#  Run it as follows: 
#     $G4BIN/Linux-g++/test22 test22.it10-verbose.in | ./process-GammaN.awk
#
# First version: John Apostolakis, 7 June 2012
#
BEGIN               { started=0;  verbose=1;     ## Not started
                      banner(); 
                    }

/PhotonInelastic/   { numPhotonA++; 
                       if( verbose > 6 ) { print $10;  }  
                     } 
/ElectroNuclear/   { numElectronA++; } 
/PositronNuclear/   { numPositronA++; } 

/gun/ || /run/ ||  /Material/  { if( started ) { report_last_energy(); clear_sums();  started=0; }  }

/beamOn/            { Events=$2 ;  started=1; } 

/energy/            { Energy= $2; Unit=$3; }

/Material/          { Material= $2; }

/gun\/particle/     { Primary= $2; } 

END                 { report_last_energy(); } 

func  banner(){
  printf  "%10s  %10s  %7s  %5s  %8s :  %8s %8s %8s \n",   "Material", "Primary",  "Energy" , "Unit",  "Events", "PhotonA", "ElectronA", "PositronA" ; 
}

func  report_last_energy(){
  printf  "%10s  %10s  %7g  %5s  %8d :  %8d %8d %8d \n",   Material, Primary,  Energy , Unit,  Events, numPhotonA, numElectronA, numPositronA; 
}

func  clear_sums() {
  numPhotonA=0 ; 
  numElectronA=0 ; 
  numPositronA=0 
}

# /run/verbose 2
# /run/initialize
 # Bias factor for Phot Nuclear is set to 100
 # Bias factor for ElectronNuclear is set to 1000
 # Bias factor for PositronNuclear is set to 1000 - same process instance is used as for electron 
# /gun/direction 1 0 0
# /tracking/verbose 1
# /mydet/SelectMaterial Pb
# /gun/particle gamma
# 
# /gun/energy 1.5 GeV
# /run/beamOn 10
#      ElectroNuclear  Models:          CHIPSElectroNuclear: Emin(GeV)=    0  Emax(GeV)= 30000
#      ElectroNuclear  Crs sctns:          ElectroNuclearXS: Emin(GeV)=    0  Emax(GeV)= 100000
#      PhotonInelastic  Models:            CHIPSGammaNuclear: Emin(GeV)=    0  Emax(GeV)= 3.5
#      PhotonInelastic  Crs sctns:            PhotoNuclearXS: Emin(GeV)=    0  Emax(GeV)= 100000
# Event number 0
#      1      199        0        0         0        0      199       199       SPhys PhotonInelastic
