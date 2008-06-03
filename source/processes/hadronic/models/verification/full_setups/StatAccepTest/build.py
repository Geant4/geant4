#!/usr/bin/python

#----------------------------------------------------------------
# This Python script has the following input parameters:
#
#   1) the Geant4 reference tag ; example:  6.2.p02
#   2) the Geant4 Physics List ;  example:  QGSP
#   3) the Calorimeter type ;     example:  PbLAr
#   4) the Particle type ;        example:  p
#   5) the beam Energy ;          example:  20GeV
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
#  iii) mainStatAccepTest.cc
#                             which is the main program.
#
# This program is called by the shell script : mainScript.sh
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

print '  REFERENCE   = ', REFERENCE
print '  PHYSICS     = ', PHYSICS
print '  CALORIMETER = ', CALORIMETER
print '  PARTICLE    = ', PARTICLE
print '  ENERGY      = ', ENERGY
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
setupFile.write( "export G4ANALYSIS_USE=1 \n" )

setupFile.write( "export CERNLIB_DIR=$DIR_INSTALLATIONS/dirCERNLIB \n" )
setupFile.write( "export PATH=$DIR_INSTALLATIONS/diriAIDA/bin:${PATH} \n" )
setupFile.write( "eval `aida-config --runtime sh` \n" )

setupFile.write( "export DIR_STAT=$DIR_INSTALLATIONS/dirStatisticalToolkit \n" )
setupFile.write( "export LD_LIBRARY_PATH=$DIR_STAT/lib:$LD_LIBRARY_PATH \n" )

setupFile.write( "export GSL_DIR=$DIR_INSTALLATIONS/dirGSL \n" )
setupFile.write( "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GSL_DIR/lib \n" )

setupFile.write( "export PATH=$PATH:$DIR_INSTALLATIONS/dirPAW \n" )

setupFile.close()

# ----------------- Write mainStatAccepTest.cc file -------------------

mainProgram = open( "mainStatAccepTest.cc", "w" )

mainProgram.write( "#include \"G4RunManager.hh\" \n" )
mainProgram.write( "#include \"G4UImanager.hh\" \n" )
mainProgram.write( "#include \"StatAccepTestDetectorConstruction.hh\" \n" )
mainProgram.write( "#include \"LHEP.hh\" \n" )
mainProgram.write( "#include \"QGSP.hh\" \n" )
mainProgram.write( "#include \"QGSP_BERT.hh\" \n" )
###mainProgram.write( "#include \"QGSP_BERT_HP.hh\" \n" )
###mainProgram.write( "#include \"QGSP_BIC_HP.hh\" \n" )
mainProgram.write( "#include \"QGSP_BIC.hh\" \n" )
mainProgram.write( "#include \"QGSC.hh\" \n" )
mainProgram.write( "#include \"QGSP_EMV.hh\" \n" )
mainProgram.write( "#include \"FTFP.hh\" \n" )
mainProgram.write( "#include \"FTFC.hh\" \n" )
###mainProgram.write( "#include \"QGS_BIC.hh\" \n" )
###mainProgram.write( "#include \"FTF_BIC.hh\" \n" )
mainProgram.write( "#include \"StatAccepTestPrimaryGeneratorAction.hh\" \n" )
mainProgram.write( "#include \"StatAccepTestEventAction.hh\" \n" )
mainProgram.write( "#include \"StatAccepTestRunAction.hh\" \n" )
mainProgram.write( "#include \"StatAccepTestTrackingAction.hh\" \n" )
mainProgram.write( "#include \"StatAccepTestStackingAction.hh\" \n" )
mainProgram.write( "#include \"StatAccepTestAnalysis.hh\" \n" )
mainProgram.write( "#include \"G4UIterminal.hh\" \n" )
mainProgram.write( "#ifdef G4UI_USE_TCSH \n" )
mainProgram.write( "  #include \"G4UItcsh.hh\" \n" )
mainProgram.write( "#endif \n" )
mainProgram.write( "#ifdef G4VIS_USE \n" )
mainProgram.write( "  #include \"StatAccepTestVisManager.hh\" \n" )
mainProgram.write( "#endif \n" )
mainProgram.write( "#ifdef G4UI_USE_XM \n" )
mainProgram.write( "#include \"G4UIXm.hh\" \n" )
mainProgram.write( "#endif \n" )
mainProgram.write( "#include \"CLHEP/Random/RanluxEngine.h\" \n" )
mainProgram.write( "int main(int argc,char** argv) { \n" )
mainProgram.write( "  CLHEP::RanluxEngine defaultEngine( 1234567, 4 ); \n" ) 
mainProgram.write( "  CLHEP::HepRandom::setTheEngine( &defaultEngine ); \n" )
mainProgram.write( "  G4int seed = time( NULL ); \n" )
mainProgram.write( "  CLHEP::HepRandom::setTheSeed( seed ); \n" )
mainProgram.write( "  G4cout << G4endl \n" )
mainProgram.write( "         << \" ===================================================== \" << G4endl \n" )
mainProgram.write( "         << \" Initial seed = \" << seed << G4endl \n" )
mainProgram.write( "	 << \" ===================================================== \" << G4endl \n" ) 
mainProgram.write( "	 << G4endl; \n" )
mainProgram.write( "  G4RunManager* runManager = new G4RunManager; \n" )
mainProgram.write( "  runManager->SetUserInitialization( new StatAccepTestDetectorConstruction ); \n" )
mainProgram.write( "  " + PHYSICS + "  *thePL = new " + PHYSICS + "; \n" )
mainProgram.write( "  //thePL->SetDefaultCutValue( 1.0*cm ); \n" ) #***LOOKHERE***
mainProgram.write( "  runManager->SetUserInitialization( thePL ); \n" )

mainProgram.write( "  runManager->SetUserAction( new StatAccepTestPrimaryGeneratorAction ); \n" )
mainProgram.write( "  runManager->SetUserAction( new StatAccepTestRunAction ); \n" )  
mainProgram.write( "  runManager->SetUserAction( new StatAccepTestEventAction ); \n" )
mainProgram.write( "  runManager->SetUserAction( new StatAccepTestTrackingAction ); \n" )

# If the beam energy is below a given threshold (in GeV) then
# no biasing is used (i.e. the StackingAction is not used).
if ( float( EnergyValue.split()[0] ) < energyThresholdInGeVNoBiasBelow ) :
    mainProgram.write( "  //runManager->SetUserAction( new StatAccepTestStackingAction ); \n" )
else :
    mainProgram.write( "  runManager->SetUserAction( new StatAccepTestStackingAction ); \n" )

#***LOOKHERE*** Switch off the histograms (leave only the ntuple).
mainProgram.write( "  StatAccepTestAnalysis::getInstance()->setIsHistogramOn( false ); \n" )

mainProgram.write( "#ifdef G4VIS_USE \n" )
mainProgram.write( "  StatAccepTestVisManager *visManager = new StatAccepTestVisManager; \n" )
mainProgram.write( "  visManager->Initialize(); \n" )
mainProgram.write( "#endif \n" )        
mainProgram.write( "  runManager->Initialize(); \n" )
mainProgram.write( "  G4UImanager* UI = G4UImanager::GetUIpointer(); \n" )   
mainProgram.write( "  if ( argc==1 ) {   // Define UI session for interactive mode. \n" )
mainProgram.write( "    G4UIsession* session = 0; \n" )
mainProgram.write( "#ifdef G4UI_USE_XM \n" )
mainProgram.write( "    session = new G4UIXm(argc,argv); \n" )
mainProgram.write( "#else \n" )
mainProgram.write( "#ifdef G4UI_USE_TCSH \n" )
mainProgram.write( "    session = new G4UIterminal(new G4UItcsh); \n" )      
mainProgram.write( "#else \n" )
mainProgram.write( "    session = new G4UIterminal(); \n" )
mainProgram.write( "#endif \n" )
mainProgram.write( "#endif \n" )
mainProgram.write( "#ifdef G4VIS_USE \n" )
mainProgram.write( "    // Create empty scene \n" )
mainProgram.write( "    G4String visCommand = \"/vis/scene/create\"; \n" )
mainProgram.write( "    UI->ApplyCommand(visCommand); \n" )
mainProgram.write( "    // Choose one default viewer (you can always change it later on) \n" ) 
mainProgram.write( "#ifdef WIN32 \n" )
mainProgram.write( "    visCommand = \"/vis/open VRML2FILE\"; \n" )
mainProgram.write( "#else \n" )
mainProgram.write( "    // visCommand = \"/vis/open VRML2\"; \n" )
mainProgram.write( "    visCommand = \"/vis/open OGLIX\"; \n" )
mainProgram.write( "#endif \n" )
mainProgram.write( "    UI->ApplyCommand(visCommand); \n" )
mainProgram.write( "    visCommand = \"/vis/viewer/flush\"; \n" )
mainProgram.write( "    UI->ApplyCommand(visCommand); \n" )
mainProgram.write( "    visCommand = \"/tracking/storeTrajectory 1\"; \n" )
mainProgram.write( "    UI->ApplyCommand(visCommand); \n" )
mainProgram.write( "#endif \n" )
mainProgram.write( "#ifdef G4UI_USE_XM \n" )
mainProgram.write( "    // Customize the G4UIXm menubar with a macro file : \n" )
mainProgram.write( "    UI->ApplyCommand(\"/control/execute gui.g4\"); \n" )
mainProgram.write( "#else \n" )
mainProgram.write( "    G4cout << \"Now, please, apply beamOn command...\" << G4endl; \n" )
mainProgram.write( "#endif \n" )
mainProgram.write( "    session->SessionStart(); \n" )
mainProgram.write( "    delete session; \n" )
mainProgram.write( "  } else {   // Batch mode \n" )
mainProgram.write( "    G4String command = \"/control/execute \"; \n" )
mainProgram.write( "    G4String fileName = argv[1]; \n" )
mainProgram.write( "    UI->ApplyCommand(command+fileName); \n" )
mainProgram.write( "  } \n" )
mainProgram.write( "  G4cout << G4endl \n" ) 
mainProgram.write( "	 << \" ===================================================== \" << G4endl \n" )
mainProgram.write( "         << \" Final random number = \" " )
mainProgram.write( "         << CLHEP::HepRandom::getTheEngine()->flat() << G4endl \n" )
mainProgram.write( "	 << \" ===================================================== \" << G4endl \n" ) 
mainProgram.write( "         << G4endl; \n" )
mainProgram.write( "  // job termination \n" )
mainProgram.write( "#ifdef G4VIS_USE \n" )
mainProgram.write( "  delete visManager; \n" )
mainProgram.write( "#endif \n" )
mainProgram.write( "  delete runManager; \n" )
mainProgram.write( "  return 0; \n" )
mainProgram.write( "} \n" )

mainProgram.close()

# ==============================================================

print '  ========== END build.py ========== '
