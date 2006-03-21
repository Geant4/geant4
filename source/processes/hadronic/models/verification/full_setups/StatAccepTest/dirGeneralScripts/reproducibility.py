#!/usr/bin/python

#----------------------------------------------------------------
# Last update: 21-Mar-2006
#
# This python script, which has no input parameters, makes 3
# tests for checking the reproducibility of the sequence of
# random numbers for a Geant4 application.
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
# are specified in ***LOOKHERE*** .
#
# The 3 tests for reproducibility of the sequence of random
# numbers, which use the same starting seed, and then look
# at the random number which is printed at the end of the
# execution of the application, are the following:
#   1) run twice  2*number_of_events ;
#   2) run first  number_of_events  and then, in the same
#      job, run another  number_of_events ;
#   3) run first  number_of_events  ; then, in another job,
#      run another  number_of_events+1  starting with the
#      seed that has been saved at the beginning of the last
#      event (currentEvent.rndm).
# The reproducibility is guaranteed if the same final random
# number is produced in each case.
#
# This script assumes that in the directory where it is run
# the following files exists:
#   -  setup.sh : Bash-shell script, which defines the
#                 environment necessary to run the executable:
#                 $G4BIN/$G4SYSTEM/mainStatAccepTest
#   -  mainStatAccepTest.cc : the main program, where the
#                             Physics List is selected;
#   -  start.rndm : the starting seed (it does not matter
#                   how you get it: you can copy from any
#                   currentEvent.rndm or currentRun.rndm
#                   files).
#
# The main result of the script is printed out on the screen.
# As by product, the script writes 3 Geant4 macro files:
# reproducibility1.g4, reproducibility2.g4, reproducibility3.g4
# and and 4 log files: out1.log, out1b.log, out2.log, out3.log;
# you do not need to look at them, and you can delete them.
# 
#----------------------------------------------------------------

import os
import sys
import string

print '  ========== START reproducibility.py ========== '

#***LOOKHERE***
ParticleType = "pi-"
###ParticleType = "pi+"
###ParticleType = "kaon-"
###ParticleType = "kaon+"
###ParticleType = "kaon0L"
###ParticleType = "neutron"
###ParticleType = "proton"

EnergyValue = "100 GeV"

Absorber = "Iron"
###Absorber = "Copper"
###Absorber = "Lead"
###Absorber = "Tungsten"
###Absorber = "PbWO4"

Active = "Scintillator"
###Active = "LiquidArgon"
###Active = "PbWO4"

isHomogeneous = "0"
###isHomogeneous = "1"

BfieldValue = "0 tesla"

NumEvents = "10"

#***endLOOKHERE***

print '  ParticleType  = ', ParticleType                 
print '  EnergyValue   = ', EnergyValue
print '  Absorber      = ', Absorber
print '  Active        = ', Active
print '  isHomogeneous = ', isHomogeneous
print '  BfieldValue   = ', BfieldValue
print '  NumEvents     = ', NumEvents

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
    g4file.write( "/gun/particle " + ParticleType + " \n" )
    g4file.write( "/gun/energy " + EnergyValue + " \n" )			
    g4file.write( "#/mydet/setField " + BfieldValue + " \n" )    
    g4file.write( "/mydet/absorberMaterial " + Absorber + " \n" )		
    g4file.write( "/mydet/activeMaterial " + Active + " \n" )
    g4file.write( "/mydet/isCalHomogeneous " + isHomogeneous + " \n" )		
    g4file.write( "/mydet/isUnitInLambda 1 \n" )		
    g4file.write( "/mydet/absorberTotalLength 10.0 \n" )	
    g4file.write( "/mydet/activeLayerNumber 100 \n" )		
    g4file.write( "/mydet/readoutLayerNumber 20 \n" )		
    g4file.write( "/mydet/activeLayerSize 4.0 \n" )		
    g4file.write( "/mydet/isRadiusUnitInLambda 1 \n" )	
    g4file.write( "/mydet/radiusBinSize 0.25 \n" )		
    g4file.write( "/mydet/radiusBinNumber 10 \n" )	
    g4file.write( "/mydet/update \n" )
    if ( i == 0 ) :
        N = str( 2 * int( NumEvents ) )
        g4file.write( "/run/beamOn " + N + " \n" )
    elif ( i == 1 ) :
        g4file.write( "/run/beamOn " + NumEvents + " \n" )
        g4file.write( "/run/beamOn " + NumEvents + " \n" )
    else :
        g4file.write( "/run/beamOn " + NumEvents + " \n" )
        g4file.write( "/random/resetEngineFrom currentEvent.rndm \n" )
        N = str( int( NumEvents ) + 1 )
        g4file.write( "/run/beamOn " + N + " \n" )
    g4file.close()

# --- Run the tests and get the final random numbers from the log files ---
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
    os.system( " . setup.sh ; mainStatAccepTest " +
               nameMacro + " > " + nameOut + " 2>&1" )
    logfile = open( nameOut, "r" )
    for line in logfile :
        if ( line.find( "Final random number" ) > - 1 ) :
            r = float( line.split("=")[1] )
            print ' ', i, ')  Final random number = ', r
            vecR[ i ] = r
    logfile.close()

# --- Analyse the results of the tests ---
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
    print '  *** REPRODUCIBILITY ?  NO  *** '
print '  ****************************** \n'

print '  ========== END reproducibility.py ========== '
