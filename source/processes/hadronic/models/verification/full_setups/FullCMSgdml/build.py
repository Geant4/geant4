#!/usr/bin/python

#----------------------------------------------------------------
# This Python script has the following input parameters:
#
#   1) the Geant4 reference tag ; example:  9.1.p02
#   2) the Geant4 Physics List ;  example:  QGSP
#   6) the Number of Events ;     example:  5k
#   7) the B field ;              example:  4T
#
# This script produces in output two files:
#
#    i) run.g4
#               which is the Geant4 command files corresponding
#               to the choice specified in the input parameters.
#   ii) setup.sh
#               which is the setup script that sets the needed
#               environmental variables.
#
# This program is called by the shell script : mainScript.sh
#----------------------------------------------------------------

import os
import sys
import string

print '  ========== START build.py ========== '

# -------------- Get input parameters and print them ------------
REFERENCE    = sys.argv[1]
PHYSICS      = sys.argv[2]
EVENTS       = sys.argv[3]
BFIELD       = sys.argv[4]

print '  REFERENCE   = ', REFERENCE
print '  PHYSICS     = ', PHYSICS
print '  EVENTS      = ', EVENTS
print '  BFIELD      = ', BFIELD

# ---------------- Release ---------------------
Release = "dirGeant4-" + REFERENCE
print '  Release = ', Release                 
                 
# ---------------- Physics ---------------------
if ( PHYSICS != "LHEP"          and
     PHYSICS != "LHEP_GN"       and
     PHYSICS != "QGSP"          and
     PHYSICS != "QGSP_GN"       and
     PHYSICS != "QGSP_BERT"     and
     PHYSICS != "QGSP_BIC"      and
     PHYSICS != "QGSP_BERT_HP"  and
     PHYSICS != "QGSP_BIC_HP"   and
     PHYSICS != "QGSC"          and
     PHYSICS != "QGSC_BERT"     and
     PHYSICS != "QGSP_EMV"      and
     PHYSICS != "FTFP"          and 
     PHYSICS != "FTFC"          and
     PHYSICS != "FTFP_BERT"     and 
     PHYSICS != "FTFP_BIC"      and
     PHYSICS != "QGS_BIC"       and
     PHYSICS != "FTF_BIC"
     
   ) :
    print '  ***ERROR*** in build.py : WRONG PHYSICS LIST = ', PHYSICS
    sys.exit( 1 )        

# ---------------- Number of events ------------
NumEvents = ""
for character in EVENTS :
    if ( character.isdigit() ) :
        NumEvents = NumEvents + character
    elif ( character == "k" ) :
        NumEvents = NumEvents + "000"
    else :
        print '  ***ERROR*** in build.py : WRONG NUMBER OF EVENTS = ', EVENTS
        sys.exit( 6 )        
print '  NumEvents = ', NumEvents

# ---------------- B field -----------------
BfieldValue = ""
isNumericPart = 1
for character in BFIELD :
    if ( isNumericPart ) :
        if ( character.isdigit()  or  character == "." ) :
            BfieldValue = BfieldValue + character
        elif ( character.isalpha() ) :
            BfieldValue = BfieldValue + " " + character
            isNumericPart = 0
    else :
        if ( character.isalpha() ) :
            BfieldValue = BfieldValue + character
        elif ( character.isdigit() ) :
            print '  ***ERROR*** in build.py : WRONG BFIELD = ', BFIELD
            sys.exit( 7 )
if ( isNumericPart ) :   # If the unit is not specified, Tesla is assumed.
    BfieldValue =  BfieldValue + " tesla" 
print '  BfieldValue = ', BfieldValue

# ==============================================================

# ----------------- Write Geant4 command file -------------------

g4file = open( "run.g4", "w" )

g4file.write( "/random/setSavingFlag 1 \n" )		
g4file.write( "/run/verbose 1 \n" )			
g4file.write( "/event/verbose 0 \n" ) 			
g4file.write( "/tracking/verbose 0 \n" )

if ( BfieldValue != ""     and
     BfieldValue != "0"    and  BfieldValue != "0 tesla"    and
     BfieldValue != "0."   and  BfieldValue != "0. tesla"   and
     BfieldValue != "0.0"  and  BfieldValue != "0.0 tesla"  ) :
    g4file.write( "/mydet/setField " + BfieldValue + " \n" )    

g4file.write( "/run/beamOn " + NumEvents + " \n" )

g4file.close()

# ----------------- Write setup file -------------------

setupFile = open( "setup.sh", "w" )

###setupFile.write( "export VO_GEANT4_SW_DIR=/users/ribon/dirGrid/dirJun08 \n" )   #***LOOKHERE***

# In the American sites, the environmental variable  $VO_GEANT4_SW_DIR
# is not defined. Its equivalent is:  $OSG_APP/geant4 .
setupFile.write( "if [ -d $VO_GEANT4_SW_DIR/dirInstallations ] ; then \n" )
setupFile.write( "  echo VO_GEANT4_SW_DIR=$VO_GEANT4_SW_DIR \n" )
setupFile.write( "else \n")
setupFile.write( "  if [ -d $OSG_APP/geant4 ] ; then \n" )
setupFile.write( "    echo OSG_APP=$OSG_APP \n" )
setupFile.write( "    echo set VO_GEANT4_SW_DIR=$OSG_APP/geant4 \n" )
setupFile.write( "    export VO_GEANT4_SW_DIR=$OSG_APP/geant4 \n" )
setupFile.write( "  else \n")
setupFile.write( "    echo ***ERROR*** : VO_GEANT4_SW_DIR and OSG_APP are undefined or unaccessible! \n" )
setupFile.write( "    exit 1 \n")
setupFile.write( "  fi \n")
setupFile.write( "fi \n")

setupFile.write( "export DIR_INSTALLATIONS=$VO_GEANT4_SW_DIR/dirInstallations \n" )

setupFile.write( "export " + PHYSICS + "=1 \n" )

setupFile.write( "export G4LIB_BUILD_GDML=1 \n" )
setupFile.write( "export G4LIB_USE_GDML=1 \n" )
setupFile.write( "export XERCESCROOT=$DIR_INSTALLATIONS/dirXERCESC \n" )
setupFile.write( "export LD_LIBRARY_PATH=$XERCESCROOT/lib:$LD_LIBRARY_PATH \n" )

###setupFile.write( "export PATH=$DIR_INSTALLATIONS/dirGCC/bin:$PATH \n" )
###setupFile.write( "export LD_LIBRARY_PATH=$DIR_INSTALLATIONS/dirGCC/lib:$LD_LIBRARY_PATH \n" )

setupFile.write( "export G4SYSTEM=Linux-g++ \n" )
setupFile.write( "export G4_RELEASE=" + Release + " \n" )

# Look for the Geant4 release in the parent directory of the current
# directory (i.e. StatAccepTest/.. ): if you find it there, use this,
# otherwise point to the installation directory.
isLocalGeant4 = 0
parentDir = os.getcwd() + "/.."
#print ' parentDir = ', parentDir
for iFile in os.listdir( parentDir ) :
    if ( iFile == Release  ) :
        isLocalGeant4 = 1
        #print ' FOUND ', Release
        break
if isLocalGeant4 :
    setupFile.write( "export G4INSTALL=" + parentDir + "/$G4_RELEASE \n" )
    setupFile.write( "export G4LIB=" + parentDir + "/$G4_RELEASE/lib \n" )
else :
    setupFile.write( "export G4INSTALL=$DIR_INSTALLATIONS/$G4_RELEASE \n" )
    setupFile.write( "export G4LIB=$DIR_INSTALLATIONS/$G4_RELEASE/lib \n" )

setupFile.write( "if [ -d $DIR_INSTALLATIONS/$G4_RELEASE/data ] ; then \n" )
setupFile.write( "  export G4LEVELGAMMADATA=$G4INSTALL/data/PhotonEvaporation \n" )
setupFile.write( "  export G4RADIOACTIVEDATA=$G4INSTALL/data/RadioactiveDecay \n" )
setupFile.write( "  export G4LEDATA=$G4INSTALL/data/G4EMLOW \n" )
# Geant4 versions >= 9.0 use G4NEUTRONHPDATA
setupFile.write( "  export G4NEUTRONHPDATA=$G4INSTALL/data/G4NDL \n")
setupFile.write( "  export NeutronHPCrossSections=$G4INSTALL/data/G4NDL \n")
setupFile.write( "else \n")
setupFile.write( "  export G4LEVELGAMMADATA=$DIR_INSTALLATIONS/dirG4DATA/PhotonEvaporation2.0 \n" )
setupFile.write( "  export G4RADIOACTIVEDATA=$DIR_INSTALLATIONS/dirG4DATA/RadioactiveDecay3.2 \n" )
setupFile.write( "  export G4LEDATA=$DIR_INSTALLATIONS/dirG4DATA/G4EMLOW5.1 \n" )
setupFile.write( "  export G4NEUTRONHPDATA=$DIR_INSTALLATIONS/dirG4DATA/G4NDL3.12 \n")
setupFile.write( "  export NeutronHPCrossSections=$DIR_INSTALLATIONS/dirG4DATA/G4NDL3.12 \n")
setupFile.write( "fi \n")
#
setupFile.write( "echo --- Data libraries --- \n")
setupFile.write( "echo G4LEVELGAMMADATA=$G4LEVELGAMMADATA \n")
setupFile.write( "echo G4RADIOACTIVEDATA=$G4RADIOACTIVEDATA \n")
setupFile.write( "echo G4LEDATA=$G4LEDATA \n")
setupFile.write( "echo NeutronHPCrossSections=$NeutronHPCrossSections \n")
setupFile.write( "echo G4NEUTRONHPDATA=$G4NEUTRONHPDATA \n")
setupFile.write( "echo ---------------------- \n")
#
if ( REFERENCE.find( "9.2") >= 0 ) :
    setupFile.write( "export CLHEP_BASE_DIR=$DIR_INSTALLATIONS/dirCLHEP-2.0.3.3 \n" )
else :
    setupFile.write( "export CLHEP_BASE_DIR=$DIR_INSTALLATIONS/dirCLHEP \n" )
#
setupFile.write( "export CLHEP_INCLUDE_DIR=$CLHEP_BASE_DIR/include \n" )
setupFile.write( "export CLHEP_LIB_DIR=$CLHEP_BASE_DIR/lib \n" )
setupFile.write( "export CLHEP_LIB=CLHEP \n" )

setupFile.write( "export G4UI_USE_TCSH=1 \n" )
setupFile.write( "export G4VIS_NONE=1 \n" )

setupFile.write( "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4LIB/$G4SYSTEM \n" )
setupFile.write( "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CLHEP_LIB_DIR \n" )
setupFile.write( "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4LIB/plists/$G4SYSTEM \n" )

setupFile.write( "export G4WORKDIR=$PWD \n" )
setupFile.write( "export PATH=$PATH:$G4WORKDIR/bin/$G4SYSTEM \n" )
###setupFile.write( "export G4ANALYSIS_USE=1 \n" )

setupFile.close()

# ==============================================================

print '  ========== END build.py ========== '
