#!/usr/bin/python

#-------------------------------------------------------------------------
# Last update: 14-Jun-2013
#
# This python script, which has at most one input parameter
# (the configuration in which Geant4 has been built: this is needed only
#  in Windows, because the executable is placed in a subdirectory whose
#  name is the configuration), tests the reproducibility of the sequence
# of random numbers for the FTFP_BERT physics list, using a Geant4
# application.
# The application consists of shooting a particle on a simple
# sampling calorimeter, for a few specified combinations of:
#   beam particle, beam energy, calorimeter type
#
# The script prints out the result on the screen.
# As by product, reproducibility.py  writes,
# for each considered case:
#   o  K Geant4 macro files:
#                   reproducibility_0.g4-xxx 
#                   reproducibility_1.g4-xxx
#                   ...
#                   reproducibility_K.g4-xxx
#   o  K log files: 
#                   out_0.log-xxx
#                   out_1.log-xxx
#                   ...
#                   out_K.log-xxx
#   o  K+1 random generator status files:
#                   currentRun.rndm    (you can delete it)
#                   currentEvent.rndm  (you can delete it)
#                   event_0.rndm-xxx
#                   event_1.rndm-xxx
#                   ...
#                   event_(K-1).rndm-xxx
# where  K = NumSingleEventChecks 
# and    "xxx"  is the case considered, i.e.
#        "particle-absorber-active"
#
# This script assumes that the environmental variables needed
# to run the application have been already properly defined;
# in particular, the executable:
#    $G4BIN/$G4SYSTEM/test68
# should exist.
# Furthermore, the script assumes that in the same directory
# where the script is run the following file exists:
#    start.rndm
# which is the starting seed.
#
#----------------------------------------------------------------

import os
import sys
import string
import subprocess
import platform

print '  ========== START reproducibility_k0lFeSci.py ========== '

# beamParticle_type: pi-, pi+, kaon-, kaon+, kaon0L, neutron, proton,
#                    anti_proton, anti_neutron, deuteron, triton, alpha
# absorber_material: Iron, Copper, Lead, Tungsten, PbWO4
# active_material:   Scintillator, LiquidArgon, PbWO4

dictionary_1stCase = { 
                      'beamParticle_type': 'kaon0L',
	              'beamParticle_kineticEnergy': '20 GeV',
                      'absorber_material': 'Iron',
                      'active_material': 'Scintillator'
                     }

listCases = [ 
              dictionary_1stCase,
#              dictionary_2ndCase,
#              dictionary_3rdCase,
#              dictionary_4thCase,
#              dictionary_5thCase,
#              dictionary_6thCase,
#              dictionary_7thCase,
#              dictionary_8thCase,
#              dictionary_9thCase,
#              dictionary_10thCase,
            ] 

BfieldValue = "4 tesla"
#
# Choose the number of events for the first loop, 
# the number of single event checks, and
# the gap between successive extra single event checks.
NumEvents = "50"
NumSingleEventChecks = 100;
GapBetweenExtraSingleEventChecks = 0;

if ( NumEvents < 0 ) :
    print ' Warning: NumEvents = ', NumEvents, '  < 0 : set it to 0 !'
    NumEvents = 0
if ( NumSingleEventChecks < 0 ) :
    print ' Warning: NumSingleEventChecks = ', NumSingleEventChecks, '  < 0 : set it to 0 !'
    NumSingleEventChecks = 0
if ( GapBetweenExtraSingleEventChecks < 0 ) :
    print ' Warning: GapBetweenExtraSingleEventChecks = ', GapBetweenExtraSingleEventChecks, '  < 0 : set it to 0 !'
    GapBetweenExtraSingleEventChecks = 0

# 10-Oct-2012 : In Windows the Geant4 configuration is used as subdirectory
#               name where to put the executable
configuration = ""
if ( len( sys.argv ) > 1 ) :
  configuration = sys.argv[1]
  print ' configuration=', configuration

# Windows platform needs a special treatment.
executable_path = "."
executable_name = "test68"
copy_command = "cp"
rndmBasename = "/currentEvent.rndm"
if ( platform.system() == 'Windows' ) :
  executable_path = configuration
  executable_name += ".exe"
  copy_command = "copy"
  rndmBasename = "\currentEvent.rndm"
  print ' executable_path=', executable_path, ' ; executable_name=', executable_name, ' ; copy_command=', copy_command, ' ; rndmBasename=', rndmBasename 

for iCase in listCases :

    iParticle = iCase[ 'beamParticle_type' ]
    iEnergy = iCase[ 'beamParticle_kineticEnergy' ]
    iAbsorber = iCase[ 'absorber_material' ]
    iActive = iCase[ 'active_material' ]

    if ( iAbsorber !=  iActive ) :
        isHomogeneous = "0"
    else :
        isHomogeneous = "1"

    print '  ------------------------------------------------------------- '
    print '  ParticleType                     = ', iParticle                 
    print '  EnergyValue                      = ', iEnergy
    print '  Absorber                         = ', iAbsorber
    print '  Active                           = ', iActive
    print '  isHomogeneous                    = ', isHomogeneous
    print '  BfieldValue                      = ', BfieldValue
    print '  NumEvents                        = ', NumEvents
    print '  NumSingleEventChecks             = ', NumSingleEventChecks
    print '  GapBetweenExtraSingleEventChecks = ', GapBetweenExtraSingleEventChecks

    # --- Write Geant4 command files ---
    suffix = '-' + iParticle  + '-' + iAbsorber + '-' + iActive
    randomEventFilename = 'Dir' + suffix + rndmBasename
    if ( NumEvents > 0  and  NumSingleEventChecks > 0 ) :
        for i in range( NumSingleEventChecks+1 ) :
            g4file = open( "reproducibility_" + str( i ) + ".g4" + suffix, "w" )
            if ( i == 0 ) :
                g4file.write( "/random/resetEngineFrom start.rndm \n" )
            else :
                g4file.write( "/random/resetEngineFrom event_" + str( i-1 ) + ".rndm" + suffix + " \n" )
            g4file.write( "/random/setSavingFlag 1 \n" )	
            g4file.write( "/random/setDirectoryName Dir" + suffix + " \n")	
            g4file.write( "/run/verbose 1 \n" )			
            g4file.write( "/event/verbose 0 \n" ) 			
            g4file.write( "/tracking/verbose 0 \n" )
            g4file.write( "/gun/particle " + iParticle + " \n" )
            g4file.write( "/gun/energy " + iEnergy + " \n" )
            g4file.write( "/mydet/setField " + BfieldValue + " \n" )    
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
                g4file.write( "/run/beamOn " + NumEvents + " \n" )
                g4file.write( "/control/shell %s %s event_0.rndm%s \n" %( copy_command, randomEventFilename, suffix ) )
                for j in range( NumSingleEventChecks-1 ) :
                    g4file.write( "/run/beamOn " + str( GapBetweenExtraSingleEventChecks+1 ) + " \n" )
                    g4file.write( "/control/shell %s %s event_%d.rndm%s \n" %( copy_command, randomEventFilename, j+1, suffix ) )
            else :
                g4file.write( "/run/beamOn 1 \n" )
            g4file.close()

    # --- Run the tests and get the last "random=" numbers from the log files
    
    if ( NumSingleEventChecks > 0 ) :
        vecRuns = []
        for i in range( NumSingleEventChecks+1 ) :
            nameMacro = "reproducibility_" + str( i ) + ".g4" + suffix
            nameOut = "out_" + str( i ) + ".log" + suffix
            outfile = open( nameOut, "w" )
            ###p = subprocess.Popen( [ os.path.join( ".", "test68" ), nameMacro ], stdout=outfile, stderr=subprocess.STDOUT )
            p = subprocess.Popen( [ os.path.join( executable_path, executable_name ), nameMacro ], 
                                  stdout=outfile, stderr=subprocess.STDOUT )
            p.wait()
            outfile.close()

            # --- Look at the last "random=" number for each run
            vecRandom = []
            logfile = open( nameOut, "r" )
            ###last_r = 0.0  # To get the value as a double number (subject to rounding)
            last_r = ""   # To get the value as a hexadecimal string (no rounding)
            for line in logfile :
                if ( line.find( "random=" ) > -1 ) :
                    ###last_r = float( line[ line.find( "random=" ):].split("=")[1] )  # double number
                    last_r = line[ line.find( "random=" ):].split("=")[1].strip()   # hexadecimal string
                elif ( line.find( "Run terminated." ) > -1 ) :
                    vecRandom.append( last_r )
            vecRuns.append( vecRandom )
            logfile.close()

    # --- Analyse the results of the tests
    if ( NumSingleEventChecks > 0 ) :
        isOK = 1
        for i in range( NumSingleEventChecks+1 ) :
            if i == 0 :
                continue
            if ( vecRuns[ 0 ][ i-1 ] == vecRuns[ i ][ 0 ] ) :
                print '  test ', str( i ) ,' OK : ', vecRuns[ 0 ][ i-1 ], vecRuns[ i ][ 0 ]
            else :
                isOK = 0
                print '  test ', str( i ), ' NO : ', vecRuns[ 0 ][ i-1 ], vecRuns[ i ][ 0 ]
        print '  ****************************** ' 
        if ( isOK ) :
            print '  *** REPRODUCIBILITY ?  YES *** '
        else :
            print >> sys.stderr, 'NO REPRODUCIBILITY'
            print '  *** REPRODUCIBILITY ?  NO  *** '
        print '  ******************************'

print '  ========== END reproducibility_k0lFeSci.py ========== '
