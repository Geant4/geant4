#!/usr/bin/python

#----------------------------------------------------------------
# This Python script has the following input parameters:
#
#   1) the Geant4 reference tag ; example:  4.6.2.ref03
#   2) the Geant4 Physics List ;  example:  QGSP
#   3) the Calorimeter type ;     example:  PbLAr
#   4) the Particle type ;        example:  p
#   5) the beam Energy ;          example:  20GeV
#   6) the Number of Events ;     example:  5k
#
# This script produces in output two files:
#
#   i) run.g4
#              which is the Geant4 command files corresponding
#              to the choice specified in the input parameters.
#  ii) setup.sh
#              which is the setup script that sets the needed
#              environmental variables.
#
# This program is called by the shell script : mainScript.sh
#----------------------------------------------------------------

import os
import sys
import string

print ' ========== START writeSetup.py ========== '

# -------------- Get input parameters and print them ------------
REFERENCE    = sys.argv[1]
PHYSICS      = sys.argv[2]
CALORIMETER  = sys.argv[3]
PARTICLE     = sys.argv[4]
ENERGY       = sys.argv[5]
EVENTS       = sys.argv[6]

print 'REFERENCE   = ', REFERENCE
print 'PHYSICS     = ', PHYSICS
print 'CALORIMETER = ', CALORIMETER
print 'PARTICLE    = ', PARTICLE
print 'ENERGY      = ', ENERGY
print 'EVENTS      = ', EVENTS

# ---------------- Release ---------------------
Release = "geant" + REFERENCE
print 'Release = ', Release                 
                 
# ---------------- Particle type ---------------
ParticleType = ""
dictParticle = { 'mu-':'mu-'   , 'mu+':'mu+'  ,
                 'e-':'e-'     , 'e+':'e+'    , 'gamma':'gamma' ,
                 'pi+':'pi+'   , 'pi-':'pi-'  ,
                 'k+':'kaon+'  , 'k-':'kaon-' , 'k0L':'kaon0L' ,
                 'n':'neutron' , 'p':'proton' }
if ( dictParticle.has_key( PARTICLE ) ) :
    ParticleType = dictParticle[ PARTICLE ]
else :
    print '***ERROR*** : WRONG PARTICLE = ', PARTICLE
    sys.exit(0)        
print 'ParticleType = ', ParticleType                 
                 
# ---------------- Beam energy -----------------
EnergyValue = ""
isNumericPart = 1
for character in ENERGY :
    if ( isNumericPart ) :
        if ( character.isdigit() ) :
            EnergyValue = EnergyValue + character
        elif ( character.isalpha() ) :
            EnergyValue = EnergyValue + " " + character
            isNumericPart = 0
    else :
        if ( character.isalpha() ) :
            EnergyValue = EnergyValue + character
        elif ( character.isdigit() ) :
            print '***ERROR*** : WRONG BEAM ENERGY = ', ENERGY
            sys.exit(0)        
print 'EnergyValue = ', EnergyValue

# ---------------- Calorimeter type ------------
dictAbsorber = { 'Fe':'Iron', 'Cu':'Copper', 'Pb':'Lead', 'W':'Tungsten',
                 'PbWO4':'PbWO4' }
dictActive = { 'Sci':'Scintillator', 'LAr':'LiquidArgon', 'PbWO4':'PbWO4' }
Absorber = "PbWO4"
Active = "PbWO4"
isHomogeneous = "1"
if ( CALORIMETER != "PbWO4" ) :
    isHomogeneous = "0"
    firstPart = ""
    secondPart = ""
    isFirstPart = 1
    for character in CALORIMETER :
        if ( len( firstPart ) > 0 and character.isupper() ) :
            isFirstPart = 0
        if ( isFirstPart ) :
            firstPart = firstPart + character
        else :
            secondPart = secondPart + character
    if ( dictAbsorber.has_key( firstPart ) ) :
        Absorber = dictAbsorber[ firstPart ]
    else :
        print '***ERROR*** : WRONG ABSORBER = ', firstPart, '  ', CALORIMETER
        sys.exit(0)        
    if ( dictActive.has_key( secondPart ) ) :
        Active = dictActive[ secondPart ]
    else :
        print '***ERROR*** : WRONG ACTIVE = ', secondPart, '  ', CALORIMETER
        sys.exit(0)        
print 'Absorber = ', Absorber
print 'Active   = ', Active

# ---------------- Number of events ------------
NumEvents = ""
for character in EVENTS :
    if ( character.isdigit() ) :
        NumEvents = NumEvents + character
    elif ( character == "k" ) :
        NumEvents = NumEvents + "000"
    else :
        print '***ERROR*** : WRONG NUMBER OF EVENTS = ', EVENTS
        sys.exit(0)        
print 'NumEvents = ', NumEvents

# ==============================================================

# ----------------- Write Geant4 command file -------------------

g4file = open( "run.g4", "w" )

g4file.write( "/random/setSavingFlag 1 \n" )		
g4file.write( "/run/verbose 1 \n" )			
g4file.write( "/event/verbose 0 \n" ) 			
g4file.write( "/tracking/verbose 0 \n" )

g4file.write( "/gun/particle " + ParticleType + " \n" )

g4file.write( "/gun/energy " + EnergyValue + " \n" )			

g4file.write( "/mydet/absorberMaterial " + Absorber + " \n" )		
g4file.write( "/mydet/activeMaterial " + Active + " \n" )
g4file.write( "/mydet/isCalHomogeneous " + isHomogeneous + " \n" )		

g4file.write( "/mydet/isUnitInLambda 1 \n" )		
g4file.write( "/mydet/absorberTotalLength 10.0 \n" )	
g4file.write( "/mydet/activeLayerNumber 100 \n" )		
g4file.write( "/mydet/activeLayerSize 4.0 \n" )		
g4file.write( "/mydet/isRadiusUnitInLambda 1 \n" )	
g4file.write( "/mydet/radiusBinSize 0.1 \n" )		
g4file.write( "/mydet/radiusBinNumber 30 \n" )	
g4file.write( "/mydet/update \n" )

g4file.write( "/run/beamOn " + NumEvents + " \n" )

g4file.close()

# ----------------- Write setup file -------------------

setupFile = open( "setup.sh", "w" )

setupFile.write( ". /afs/cern.ch/sw/geant4/dev/scripts/gcc-alt.sh 3.2.3 \n" )

setupFile.write( "export G4SYSTEM=Linux-g++ \n" )

setupFile.write( "RELEASE=" + Release + " \n" )
setupFile.write( "PLATFORM=rh73_gcc323/ \n" )
setupFile.write( "DIR_SPECIFIC=/afs/cern.ch/sw/geant4/releases/specific/$PLATFORM \n" )
setupFile.write( "export G4INSTALL=/afs/cern.ch/sw/geant4/releases/share/$RELEASE \n" )
setupFile.write( "export G4LIB=$DIR_SPECIFIC/$RELEASE/lib \n" )

setupFile.write( "export G4LEVELGAMMADATA=/afs/cern.ch/sw/geant4/dev/data/PhotonEvaporation \n" )
setupFile.write( "export G4RADIOACTIVEDATA=/afs/cern.ch/sw/geant4/dev/data/RadiativeDecay \n" )
setupFile.write( "export G4LEDATA=/afs/cern.ch/sw/geant4/dev/data/G4EMLOW \n" )
setupFile.write( "export NeutronHPCrossSections=/afs/cern.ch/sw/geant4/dev/data/G4NDL \n")

setupFile.write( "export CLHEP_BASE_DIR=/afs/cern.ch/sw/geant4/dev/CLHEP/1.8.1.0/redhat73_gcc323 \n" )
setupFile.write( "export CLHEP_INCLUDE_DIR=$CLHEP_BASE_DIR/include \n" )
setupFile.write( "export CLHEP_LIB_DIR=$CLHEP_BASE_DIR/lib \n" )
setupFile.write( "export CLHEP_LIB=CLHEP \n" )

setupFile.write( "export G4ANALYSIS_USE=1 \n" )
setupFile.write( "export G4UI_USE_TCSH=1 \n" )

setupFile.write( "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4LIB/$G4SYSTEM \n" )
setupFile.write( "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4LIB/.lists_build/$G4SYSTEM \n" )
setupFile.write( "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CLHEP_LIB_DIR \n" )

setupFile.write( "export G4WORKDIR=$PWD \n" )
setupFile.write( "export PATH=$PATH:$G4WORKDIR/bin/$G4SYSTEM \n" )
setupFile.write( "export G4ANALYSIS_USE=1 \n" )

setupFile.write( "export PI_DIR=/afs/cern.ch/sw/lcg/app/releases/PI/PI_1_2_3/rh73_gcc323 \n" )
setupFile.write( "export PATH=$PI_DIR/bin:$PATH \n" )
setupFile.write( "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PI_DIR/lib \n" )
setupFile.write( "eval `aida-config --runtime sh` \n" )

setupFile.write( "export LD_LIBRARY_PATH=/afs/cern.ch/user/r/ribon/ExtraSpace/myPiDir/lib:$LD_LIBRARY_PATH \n" )

setupFile.write( "export GSL_DIR=/afs/cern.ch/sw/lcg/external/GSL/1.4/rh73_gcc323 \n" )
setupFile.write( "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GSL_DIR/lib \n" )

setupFile.close()

# ==============================================================

print ' ========== END writeSetup.py ========== '
