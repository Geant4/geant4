#!/usr/bin/python

#-------------------------------------------------------------------
# Last update: 01-Aug-2008
#
# This script has 1 input argument, and should be run as:
#
#         $  python printInfoLogfile.py file
#
# where the input  file  should be one of the log-files
# produced by running the  DetectorSlice  simulations.
# This script prints out (in the screen) the following
# observables:
#
#   o  CPU times
#   o  energy in the Tracker 
#   o  total visible energy in the EM calorimeter
#   o  total energy in the whole EM calorimeter
#   o  total visible energy in the HAD calorimeter
#   o  total energy in the whole HAD calorimeter
#   o  energy in the Muon detector
#   o  energy resolution in the EM calorimeter
#   o  energy resolution in the HAD calorimeter
#   o  visible energy in the EM calorimeter due to electrons
#   o  visible energy in the EM calorimeter due to muons
#   o  visible energy in the EM calorimeter due to pions
#   o  visible energy in the EM calorimeter due to kaons
#   o  visible energy in the EM calorimeter due to protons
#   o  visible energy in the EM calorimeter due to nuclei
#   o  visible energy in the HAD calorimeter due to electrons
#   o  visible energy in the HAD calorimeter due to muons
#   o  visible energy in the HAD calorimeter due to pions
#   o  visible energy in the HAD calorimeter due to kaons
#   o  visible energy in the HAD calorimeter due to protons
#   o  visible energy in the HAD calorimeter due to nuclei
#
# as extracted from the log-files.
# In the case that the  file  does not exist, or it exists
# but is not a log-file obtained by completing a running
# of DetectorSlice, then "0" is print for each of the above
# observables.
#
#-------------------------------------------------------------------

import os
import sys
import string
import math


#===============================================
#================= FUNCTIONS =================== 
#===============================================

def funExtract( theFile ) :
    # Given the file in input, it returns all the observables values,
    # in the form of pairs: value error .
    
    statusEndRun = 0
    if ( theFile ) :
        statusSubdetector = 0
        statusEMcalo = 0
        statusHADcalo = 0
        statusParticleTypeEMcalo = 0
        statusParticleTypeHADcalo = 0
        name = ""
        for line in theFile :
            #print line
            if ( line.find( "Run Summary" ) > -1 ) :
                statusEndRun = 1
            elif ( statusEndRun ) :                    
                if ( line.find( "User=" ) > -1 ) :
                    print "  cpuTime = ", \
                          float( line.split( "s " )[0].split( "=" )[1] )
                elif ( line.find( "deposited in the subdet" ) > -1 ) :
                    statusSubdetector = 1
                elif( line.find( "for EM Calo" ) > -1 ) :
                    statusEMcalo = 1
                    statusSubdetector = 0
                elif ( line.find( "for HAD Calo" ) > -1 ) :
                    statusHADcalo = 1
                    statusEMcalo = 0
                elif ( line.find( "main particle types [MeV] in EM Calo" ) > -1 ) :
                    statusParticleTypeEMcalo = 1
                    statusHADcalo = 0
                elif ( line.find( "main particle types [MeV] in HAD Calo" ) > -1 ) :
                    statusParticleTypeHADcalo = 1
                    statusParticleTypeEMcalo = 0

                elif ( statusSubdetector ) :
                    if ( line.find( "Tracker" ) > -1 ) :
                        print "  tracker_energy = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "EM Calo active layers" ) > -1 ) :
                        print " em_visEnergy = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "EM Calo" ) > -1  and
                           line.find( "active layers" ) == -1 ) :
                        print " em_totEnergy = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "HAD Calo active layers" ) > -1 ) :
                        print " had_visEnergy = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "HAD Calo" ) > -1  and
                           line.find( "active layers" ) == -1 ) :
                        print " had_totEnergy = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "Muon detector" ) > -1 ) :
                        print " muon_detector_energy = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )

                elif ( statusEMcalo ) :
                    if ( line.find( "energy resolution" ) > -1 ) :
                        resolutionValue = float( line.split( "=" )[1].split( "+/-" )[0] )
                        resolutionValue *= 100.0
                        print "  em_resolution = %.2f" % resolutionValue

                elif ( statusHADcalo ) :
                    if ( line.find( "energy resolution" ) > -1 ) :
                        resolutionValue = float( line.split( "=" )[1].split( "+/-" )[0] )
                        resolutionValue *= 100.0
                        print "  had_resolution = %.2f" % resolutionValue
    
                elif ( statusParticleTypeEMcalo ) :
                    if ( line.find( "electron" ) > -1 ) :
                        name = "em_electronEvis"
                    elif ( line.find( "muon" ) > -1 ) :
                        name = "em_muonEvis"
                    elif ( line.find( "pion" ) > -1 ) :
                        name = "em_pionEvis"
                    elif ( line.find( "kaon" ) > -1 ) :
                        name = "em_kaonEvis"
                    elif ( line.find( "proton" ) > -1 ) :
                        name = "em_protonEvis"
                    elif ( line.find( "nuclei" ) > -1 ) :
                        name = "em_nucleusEvis"
                    elif ( line.find( "<E_vis> =" ) > -1 ) :
                        print "  ", name, " = %.1f" % \
                              float( line.split( "(" )[1].split( "+/-" )[0] )

                elif ( statusParticleTypeHADcalo ) :
                    if ( line.find( "electron" ) > -1 ) :
                        name = "had_electronEvis"
                    elif ( line.find( "muon" ) > -1 ) :
                        name = "had_muonEvis"
                    elif ( line.find( "pion" ) > -1 ) :
                        name = "had_pionEvis"
                    elif ( line.find( "kaon" ) > -1 ) :
                        name = "had_kaonEvis"
                    elif ( line.find( "proton" ) > -1 ) :
                        name = "had_protonEvis"
                    elif ( line.find( "nuclei" ) > -1 ) :
                        name = "had_nucleusEvis"
                    elif ( line.find( "<E_vis> =" ) > -1 ) :
                        print "  ", name, " = %.1f" % \
                              float( line.split( "(" )[1].split( "+/-" )[0] )

    if ( statusEndRun == 0 ) :
        print "  cpuTime =  0  "
        print "  tracker_energy =  0  "
        print "  em_visEnergy =  0  "
        print "  em_totEnergy =  0  "
        print "  had_visEnergy =  0  "
        print "  had_totEnergy =  0  "
        print "  muon_detector_energy =  0  "
        print "  em_resolution =  0  "
        print "  em_electronEvis =  0  "
        print "  em_muonEvis =  0  "
        print "  em_pionEvis =  0  "
        print "  em_kaonEvis =  0  "
        print "  em_protonEvis =  0  "
        print "  em_nucleusEvis =  0  "
        print "  had_resolution =  0  "
        print "  had_electronEvis =  0  "
        print "  had_muonEvis =  0  "
        print "  had_pionEvis =  0  "
        print "  had_kaonEvis =  0  "
        print "  had_protonEvis =  0  "
        print "  had_nucleusEvis =  0  "

    return


#===============================================
#==================== MAIN ===================== 
#===============================================

# Check whether the number of arguments are correct or not.
if ( len( sys.argv ) != 2 ) :
    print " Usage:  printInfoLogfile.py file"
else :
    if ( os.path.isfile( sys.argv[1] ) ) :
        theFile = open( sys.argv[1], "r" )
        funExtract( theFile )
        theFile.close()
    else :
        funExtract( 0 )

#-------------------------------------------------------------------

