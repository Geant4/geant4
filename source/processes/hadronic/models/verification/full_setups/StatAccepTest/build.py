#!/usr/bin/python

#----------------------------------------------------------------
# Last update: 01-Dec-2009.
#
# This Python script has the following input parameters:
#
#   1) the Geant4 reference tag ; example:  9.2.p01
#   2) the Geant4 Physics List ;  example:  QGSP
#   3) the Calorimeter type ;     example:  PbLAr
#   4) the Particle type ;        example:  p
#   5) the beam Energy ;          example:  20GeV
#   6) the Number of Events ;     example:  5k
#   7) the B field ;              example:  4T
#   8) the string `uname -i`;     example:  i386 (or x86_64)
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
#
# 19-Mar-2009: this script does not write the main source file
#              (mainStatAccepTest.cc) anymore, because the
# executables are pre-built and distributed with the application,
# avoiding complication with the local environment.
#
#----------------------------------------------------------------

import os
import sys
import string

print '  ========== START build.py ========== '

# Below a certain beam energy threshold (in GeV) biasing is not
# applied and also the default range production cut (700 mu) is used. 
energyThresholdInGeVNoBiasBelow = 10.0   #***LOOKHERE***

print '  energy Threshold No Bias Below = ', energyThresholdInGeVNoBiasBelow, ' GeV'

# -------------- Get input parameters and print them ------------
REFERENCE    = sys.argv[1]
PHYSICS      = sys.argv[2]
CALORIMETER  = sys.argv[3]
PARTICLE     = sys.argv[4]
ENERGY       = sys.argv[5]
EVENTS       = sys.argv[6]
BFIELD       = sys.argv[7]
UNAMEI       = sys.argv[8]

print '  REFERENCE   = ', REFERENCE
print '  PHYSICS     = ', PHYSICS
print '  CALORIMETER = ', CALORIMETER
print '  PARTICLE    = ', PARTICLE
print '  ENERGY      = ', ENERGY
print '  EVENTS      = ', EVENTS
print '  BFIELD      = ', BFIELD
print '  UNAMEI      = ', UNAMEI

# ---------------- Release ---------------------
Release = "dirGeant4-" + REFERENCE
print '  Release = ', Release                 

ReleaseUnameI = Release + "_" + UNAMEI
print '  ReleaseUnameI = ', ReleaseUnameI

# ---------------- Physics ---------------------
if ( PHYSICS != "LHEP"            and
     PHYSICS != "LHEP_GN"         and
     PHYSICS != "QGSP"            and
     PHYSICS != "QGSP_GN"         and
     PHYSICS != "QGSP_EMV"        and
     PHYSICS != "QGSP_BERT"       and
     PHYSICS != "QGSP_BERT_EMV"   and
     PHYSICS != "QGSP_BIC"        and
     PHYSICS != "QGSP_BERT_HP"    and
     PHYSICS != "QGSP_BIC_HP"     and
     PHYSICS != "QGSC"            and
     PHYSICS != "QGSC_EFLOW"      and
     PHYSICS != "QGSC_BERT"       and
     PHYSICS != "QGSC_CHIPS"      and
     PHYSICS != "QGSC_QGSC"       and
     PHYSICS != "FTFP"            and 
     PHYSICS != "FTFC"            and
     PHYSICS != "FTFP_BERT"       and 
     PHYSICS != "FTFP_BIC"        and
     PHYSICS != "QGS_BIC"         and
     PHYSICS != "FTF_BIC"         and
     PHYSICS != "CHIPS"           and
     PHYSICS != "FTFP_BERT_TRV"   and
     PHYSICS != "QGSP_FTFP_BERT"  and 
     PHYSICS != "QGSP_BERT_TRV"   and 
     PHYSICS != "QGSP_INCL_ABLA"
   ) :
    print '  ***ERROR*** in build.py : WRONG PHYSICS LIST = ', PHYSICS
    sys.exit( 1 )        

# ---------------- Particle type ---------------
ParticleType = ""
dictParticle = { 'mu-':'mu-'   , 'mu+':'mu+'  ,
                 'e-':'e-'     , 'e+':'e+'    , 'gamma':'gamma' ,
                 'pi+':'pi+'   , 'pi-':'pi-'  ,
                 'k+':'kaon+'  , 'k-':'kaon-' , 'k0L':'kaon0L' ,
                 'n':'neutron' , 'p':'proton' ,
                 'nbar':'anti_neutron' , 'pbar':'anti_proton' ,
                 'd':'deuteron' , 't':'triton', 'alpha':'alpha' }
if ( dictParticle.has_key( PARTICLE ) ) :
    ParticleType = dictParticle[ PARTICLE ]
else :
    print '  ***ERROR*** in build.py : WRONG PARTICLE = ', PARTICLE
    sys.exit( 2 )        
print '  ParticleType = ', ParticleType                 
                 
# ---------------- Beam energy -----------------
EnergyValue = ""
isNumericPart = 1
for character in ENERGY :
    if ( isNumericPart ) :
        if ( character.isdigit()  or  character == "." ) :
            EnergyValue = EnergyValue + character
        elif ( character.isalpha() ) :
            EnergyValue = EnergyValue + " " + character
            isNumericPart = 0
    else :
        if ( character.isalpha() ) :
            EnergyValue = EnergyValue + character
        elif ( character.isdigit() ) :
            print '  ***ERROR*** in build.py : WRONG BEAM ENERGY = ', ENERGY
            sys.exit( 3 )        
print '  EnergyValue = ', EnergyValue

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
        print '  ***ERROR*** in build.py : WRONG ABSORBER = ', firstPart, '  ', CALORIMETER
        sys.exit( 4 )        
    if ( dictActive.has_key( secondPart ) ) :
        Active = dictActive[ secondPart ]
    else :
        print '  ***ERROR*** in build.py : WRONG ACTIVE = ', secondPart, '  ', CALORIMETER
        sys.exit( 5 )        
print '  Absorber = ', Absorber
print '  Active   = ', Active

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

###g4file.write( "/run/setCut 1.0 cm \n" ) #***LOOKHERE***

g4file.write( "/gun/particle " + ParticleType + " \n" )

g4file.write( "/gun/energy " + EnergyValue + " \n" )			

if ( BfieldValue != ""     and
     BfieldValue != "0"    and  BfieldValue != "0 tesla"    and
     BfieldValue != "0."   and  BfieldValue != "0. tesla"   and
     BfieldValue != "0.0"  and  BfieldValue != "0.0 tesla"  ) :
    g4file.write( "/mydet/setField " + BfieldValue + " \n" )    

g4file.write( "/mydet/absorberMaterial " + Absorber + " \n" )		
g4file.write( "/mydet/activeMaterial " + Active + " \n" )
g4file.write( "/mydet/isCalHomogeneous " + isHomogeneous + " \n" )		

g4file.write( "/mydet/isUnitInLambda 1 \n" )		
g4file.write( "/mydet/absorberTotalLength 10.0 \n" )	
g4file.write( "/mydet/calorimeterRadius 5.0 \n" )	
g4file.write( "/mydet/activeLayerNumber 100 \n" )		
g4file.write( "/mydet/readoutLayerNumber 20 \n" )		
g4file.write( "/mydet/activeLayerSize 4.0 \n" )		
g4file.write( "/mydet/isRadiusUnitInLambda 1 \n" )	
g4file.write( "/mydet/radiusBinSize 0.25 \n" )		
g4file.write( "/mydet/radiusBinNumber 10 \n" )	
g4file.write( "/mydet/update \n" )

g4file.write( "/run/beamOn " + NumEvents + " \n" )

g4file.close()

# ----------------- Write setup file -------------------

setupFile = open( "setup.sh", "w" )

setupFile.write( "#!/bin/sh \n" )

###setupFile.write( "export VO_GEANT4_SW_DIR=/users/ribon/dirGrid/dirDec09 \n" )   #***LOOKHERE***

# In the American sites, the environmental variable  $VO_GEANT4_SW_DIR
# is not defined. Its equivalent is:  $OSG_APP/geant4 .
setupFile.write( "if [ -d $VO_GEANT4_SW_DIR/dirInstallations ] ; then \n" )
setupFile.write( "  echo VO_GEANT4_SW_DIR=$VO_GEANT4_SW_DIR \n" )
setupFile.write( "  export DIR_INSTALLATIONS=$VO_GEANT4_SW_DIR/dirInstallations \n" )
setupFile.write( "else \n")
setupFile.write( "  if [ -d $OSG_APP/geant4 ] ; then \n" )
setupFile.write( "    echo OSG_APP=$OSG_APP \n" )
setupFile.write( "    echo set VO_GEANT4_SW_DIR=$OSG_APP/geant4 \n" )
setupFile.write( "    export VO_GEANT4_SW_DIR=$OSG_APP/geant4 \n" )
setupFile.write( "    export DIR_INSTALLATIONS=$VO_GEANT4_SW_DIR/dirInstallations_`uname -i` \n" )
setupFile.write( "  else \n")
setupFile.write( "    echo ***ERROR*** : VO_GEANT4_SW_DIR and OSG_APP are undefined or unaccessible! \n" )
setupFile.write( "    exit 1 \n")
setupFile.write( "  fi \n")
setupFile.write( "fi \n")

#***LOOKHERE***
###setupFile.write( "export PATH=$DIR_INSTALLATIONS/dirGCC/bin:$PATH \n" )
###setupFile.write( "export LD_LIBRARY_PATH=$DIR_INSTALLATIONS/dirGCC/lib:$LD_LIBRARY_PATH \n" )

setupFile.write( "export G4SYSTEM=Linux-g++ \n" )
setupFile.write( "export G4_RELEASE=" + Release + " \n" )

# Look for the Geant4 release in the parent directory of the current
# directory (i.e. StatAccepTest/.. ): if you find it there, use this,
# otherwise point to the installation directory.
isLocalGeant4 = 0
isReleaseUnameI = 0
parentDir = os.getcwd() + "/.."
#print ' parentDir = ', parentDir
for iFile in os.listdir( parentDir ) :
    if ( iFile == Release  or  iFile == ReleaseUnameI ) :
        #print ' iFile=' , iFile
        isLocalGeant4 = 1
        if ( iFile == ReleaseUnameI ) :
            isReleaseUnameI = 1
            print ' FOUND local ', ReleaseUnameI
        else :
            print ' FOUND local ', Release
        break
if isLocalGeant4 :
    if ( isReleaseUnameI ) :
        setupFile.write( "export G4INSTALL=" + parentDir + "/${G4_RELEASE}_" +
                         UNAMEI + " \n" )
        setupFile.write( "export G4LIB=" + parentDir + "/${G4_RELEASE}_" +
                         UNAMEI + "/lib \n" )
    else :
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
setupFile.write( "  export G4LEDATA=$DIR_INSTALLATIONS/dirG4DATA/G4EMLOW6.9 \n" )
setupFile.write( "  export G4NEUTRONHPDATA=$DIR_INSTALLATIONS/dirG4DATA/G4NDL3.13 \n")
setupFile.write( "  export NeutronHPCrossSections=$DIR_INSTALLATIONS/dirG4DATA/G4NDL3.13 \n")
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
if ( REFERENCE.find( "9.2" ) > -1 ) :
    setupFile.write( "export CLHEP_BASE_DIR=$DIR_INSTALLATIONS/dirCLHEP-2.0.4.2 \n" )
else :
    setupFile.write( "export CLHEP_BASE_DIR=$DIR_INSTALLATIONS/dirCLHEP-2.0.4.4 \n" )
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
setupFile.write( "export G4ANALYSIS_USEROOT=1 \n" )

setupFile.write( "export ROOTSYS=$DIR_INSTALLATIONS/dirROOT \n" )
setupFile.write( "export PATH=${PATH}:$ROOTSYS/bin \n" )
setupFile.write( "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib \n" )

# The following is needed by ROOT to find libG4TypesDict.so, which has
# been created with the reference version of Geant4 and copied by hand
# from $G4WORKDIR/tmp/$G4SYSTEM/mainStatAccepTest/ .
setupFile.write( "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4WORKDIR/lib_`uname -i`/$G4SYSTEM \n" )

# If the beam energy is below a given threshold (in GeV) then
# no biasing is used (i.e. the StackingAction is not used).
if ( float( EnergyValue.split()[0] ) < energyThresholdInGeVNoBiasBelow ) :
    setupFile.write( "unset IS_BIASING_ACTIVE \n" )
else :
    setupFile.write( "export IS_BIASING_ACTIVE=1 \n" )

#***LOOKHERE*** Switch on the extra histograms (else only the ntuple).
###setupFile.write( "export HISTOGRAMS_ON=1 \n" )

if ( REFERENCE.find( "9.2" ) > -1 ) :
    if ( PHYSICS=="CHIPS" ) :
        setupFile.write( "export PHYSLIST=QGSC_BERT \n" )
    elif ( PHYSICS=="QGSC_CHIPS" ) :
        setupFile.write( "export PHYSLIST=QGSC_BERT \n" )
    elif ( PHYSICS=="QGSC_QGSC" ) :
        setupFile.write( "export PHYSLIST=QGSC_BERT \n" )
    elif ( PHYSICS=="FTFP_BERT_TRV" ) :
        setupFile.write( "export PHYSLIST=FTFP_BERT \n" )
    elif ( PHYSICS=="QGSP_FTFP_BERT" ) :
        setupFile.write( "export PHYSLIST=QGSP_BERT \n" )
    elif ( PHYSICS=="QGSP_BERT_TRV" ) :
        setupFile.write( "export PHYSLIST=QGSP_BERT \n" )
    elif ( PHYSICS=="QGSP_INCL_ABLA" ) :
        setupFile.write( "export PHYSLIST=QGSP_BERT \n" )
    else :
        setupFile.write( "export PHYSLIST=" + PHYSICS + "\n" )
else :
    setupFile.write( "export PHYSLIST=" + PHYSICS + "\n" )

setupFile.close()

# ==============================================================

print '  ========== END build.py ========== '
