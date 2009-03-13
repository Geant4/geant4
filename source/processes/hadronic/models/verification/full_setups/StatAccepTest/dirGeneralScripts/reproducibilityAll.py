#!/usr/bin/python

#----------------------------------------------------------------
# Last update: 13-Mar-2009
#
# This python script, which has no input parameters, makes 3
# tests for checking the reproducibility of the sequence of
# random numbers for a Geant4 application, for several cases
# (beam particle, beam energy, calorimeter type).
#
# The application consists of shooting a particle on a simple
# sampling calorimeter. All its parameters:
#   - beam particle type
#   - beam particle kinetic energy
#   - absorber material
#   - active material
#   - size of the box (either in mm or in lambda of the absorber)
#   - number of layers
#   - size and number of radial rings
#   - magnetic field
#   - number of events
# are specified in ***LOOKHERE*** : for those defined as lists,
# a loop is made over each element.
#
# The 3 tests for reproducibility of the sequence of random
# numbers, which use the same starting seed, and then look
# at the random number which is printed at the end of the
# execution of the application, are the following:
#   1) run twice  number_of_events_1 + number_of_events_2 ;
#   2) run first  number_of_events_1  and then, in the same
#      job, run another  number_of_events_2 ;
#   3) run first  number_of_events_1  ; then, in another job,
#      run another  number_of_events_2 + 1  starting with the
#      seed that has been saved at the beginning of the last
#      event (currentEvent.rndm).
# The reproducibility is guaranteed if the same final random
# number is produced in each case.
#
# This script assumes that the environmental variables needed
# to run the application have been already properly defined;
# in particular, the executable:
#    $G4BIN/$G4SYSTEM/mainStatAccepTest-PHYSICS_LIST
# should exist
# (the Physics List can be selected in the main program:
#    mainStatAccepTest.cc
#  and then you should set "PHYSICS_LIST" in this script,
#  see "***LOOKHERE***" below).
# Furthermore, the script assumes that in the same directory
# where the script is run the following file exists:
#    start.rndm
# which is the starting seed
# (it does not matter how you get it: you can copy from any
#  currentEvent.rndm or currentRun.rndm files).
#
# The main result of the script is printed out on the screen.
# As by product, the script writes 3 Geant4 macro files:
# reproducibility1.g4, reproducibility2.g4, reproducibility3.g4
# and and 4 log files: out1.log, out1b.log, out2.log, out3.log;
# you do not need to look at them, and you can delete them.
# 
# NB) This script has been obtained, after few simple changes,
#     from the similar script  reproducibility.py , which
#     considers a single case instead of looping over all
#     possible cases.
#
# NB) When using this script to test reproducibility, it is
#     convenient to activate biasing, in order to test more
#     uniformily the high-energy models and the low-energy ones.
#     To do so, you have to do the following:
#       -  uncomment the following line in the main program
#          (mainStatAccepTest.cc):
#          //runManager->SetUserAction(new StatAccepTestStackingAction);
#       -  check whether the default biasing parameters
#          (see "***LOOKHERE***" in src/StatAccepTestStackingAction.cc)
#          are ok for you.
#          We suggest to change only one parameter:
#            neutronKillingEnergyThreshold( 0.0 )
#          (to avoid to kill neutrons below 1 MeV).
#     Using biasing, with the typical 100 GeV beam energy, it is
#     possible to increase by an order of magnitude the number of
#     events. A reasonable choice is:  NumEvents1 = "1000"
#                                      NumEvents2 = "100"
#
#----------------------------------------------------------------

import os
import sys
import string

print '  ========== START reproducibilityAll.py ========== '

#***LOOKHERE***

PHYSICS_LIST = 'LHEP'

listParticleTypes = [
    "pi-" ,
    "pi+" ,
    "kaon-" ,
    "kaon+" ,
    "kaon0L" , 
    "neutron" , 
    "proton" ,
    "anti_proton" , 
    "anti_neutron" ,
    "deuteron" ,
    "triton" ,
    "alpha"
    ]
#
listEnergyValues = [
    "100 GeV"
    ]
#
listAbsorbers = [
    "Iron" ,
    "Copper" , 
    "Lead" , 
    "Tungsten" , 
    "PbWO4"
    ]
#
listActiveLayer = [
    "Scintillator" ,
    "LiquidArgon"  ,
#    "PbWO4"
    ]
#
BfieldValue = "4 tesla"
#
NumEvents1 = "1000"
NumEvents2 = "100"

#***endLOOKHERE***

isOKforAllCases = 1
for iAbsorber in listAbsorbers :
    for iActive in listActiveLayer :
        if ( iAbsorber !=  iActive ) :
            isHomogeneous = "0"
        else :
            isHomogeneous = "1"
        for iParticle in listParticleTypes :
            for iEnergy in listEnergyValues :

                print '  ---------------------------------- '
                print '  ParticleType  = ', iParticle                 
                print '  EnergyValue   = ', iEnergy
                print '  Absorber      = ', iAbsorber
                print '  Active        = ', iActive
                print '  isHomogeneous = ', isHomogeneous
                print '  BfieldValue   = ', BfieldValue
                print '  NumEvents1    = ', NumEvents1
                print '  NumEvents2    = ', NumEvents2

                # --- Write Geant4 command files ---
                for i in range(3) :
                    if ( i == 0 ) :
                        g4file = open( "reproducibility1.g4", "w" )
                    elif ( i == 1 ) :
                        g4file = open( "reproducibility2.g4", "w" )
                    else :
                        g4file = open( "reproducibility3.g4", "w" )        
                    g4file.write( "/random/resetEngineFrom start.rndm \n" )
                    g4file.write( "/random/setSavingFlag 1 \n" )		
                    g4file.write( "/run/verbose 1 \n" )			
                    g4file.write( "/event/verbose 0 \n" ) 			
                    g4file.write( "/tracking/verbose 0 \n" )
                    g4file.write( "/gun/particle " + iParticle + " \n" )
                    g4file.write( "/gun/energy " + iEnergy + " \n" )
                    g4file.write( "#/mydet/setField " + BfieldValue + " \n" )    
                    g4file.write( "/mydet/absorberMaterial " + iAbsorber + " \n" )
                    g4file.write( "/mydet/activeMaterial " + iActive + " \n" )
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
                    if ( i == 0 ) :
                        N = str( int( NumEvents1 ) + int( NumEvents2 ) )
                        g4file.write( "/run/beamOn " + N + " \n" )
                    elif ( i == 1 ) :
                        g4file.write( "/run/beamOn " + NumEvents1 + " \n" )
                        g4file.write( "/run/beamOn " + NumEvents2 + " \n" )
                    else :
                        g4file.write( "/run/beamOn " + NumEvents1 + " \n" )
                        g4file.write( "/random/resetEngineFrom currentEvent.rndm \n" )
                        N = str( int( NumEvents2 ) + 1 )
                        g4file.write( "/run/beamOn " + N + " \n" )
                    g4file.close()

                # --- Run the tests and get the final random numbers from the log files
                vecR = [ 0.0 , 0.0 , 0.0 , 0.0 ]
                for i in range(4) :
                    nameMacro = "reproducibility1.g4"
                    nameOut = "out1.log"
                    if ( i == 1 ) :
                        nameOut = "out1b.log"
                    if ( i == 2 ) :
                        nameMacro = "reproducibility2.g4"
                        nameOut = "out2.log"
                    elif ( i == 3 ) :
                        nameMacro = "reproducibility3.g4"
                        nameOut = "out3.log"

                    print ' nameMacro=', nameMacro, ' nameOut=', nameOut
                    
                    os.system( " mainStatAccepTest-" + PHYSICS_LIST +
                               " " + nameMacro + " > " + nameOut + " 2>&1" )

                    logfile = open( nameOut, "r" )
                    for line in logfile :
                        if ( line.find( "Final random number" ) > - 1 ) :
                            r = float( line.split("=")[1] )
                            print ' ', i, ')  Final random number = ', r
                            vecR[ i ] = r
                    logfile.close()

                # --- Analyse the results of the tests
                isOK = 1
                if ( vecR[ 0 ] == vecR[ 1 ] ) :
                    print '  1st test OK : ', vecR[ 0 ], vecR[ 1 ]
                else :
                    isOK = 0
                    print '  1st test FAILS : ', vecR[ 0 ], vecR[ 1 ]
                if ( vecR[ 0 ] == vecR[ 2 ] ) :
                    print '  2nd test OK : ', vecR[ 0 ], vecR[ 2 ]
                else :
                    isOK = 0
                    print '  2nd test FAILS : ', vecR[ 0 ], vecR[ 2 ]
                if ( vecR[ 0 ] == vecR[ 3 ] ) :
                    print '  3rd test OK : ', vecR[0], vecR[ 3 ]
                else :
                    isOK = 0    
                    print '  3rd test FAILS : ', vecR[0], vecR[ 3 ]

                print '\n  ****************************** ' 
                if ( isOK ) :
                    print '  *** REPRODUCIBILITY ?  YES *** '
                else :
                    isOKforAllCases = 0 
                    print '  *** REPRODUCIBILITY ?  NO  *** '
                print '  ****************************** \n'

print '   '
print '   '
print ' === SUMMARY === : REPRODUCIBILITY ',
if ( isOKforAllCases ) :
    print ' OK'
else :
    print ' VIOLATED!'
print '  '

print '  ========== END reproducibilityAll.py ========== '
