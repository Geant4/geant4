#!/usr/bin/python

#-------------------------------------------------------------------
# Last update: 02-Nov-2007
#
# This script has 1 input argument, and should be run as:
#
#         $  python printInfoLogfile.py file
#
# where the input  file  should be one of the log-files
# produced by running the StatAccepTest simulations.
# This script prints out (in the screen) the following
# observables (in the case of a proper "file", otherwise
# if it is not a log-file of StatAccepTest nothing will happen):
#
#   o  CPU times
#   o  total visible energy
#   o  total energy in the whole calorimeter
#   o  visible energy in each layer (in MeV)
#   o  fractions in 4 longitudinal sections
#   o  visible energy in each ring (in MeV)
#   o  fractions in 3 radial sections
#   o  energy resolution
#   o  shower shapes per particle types
#   o  information on the number of steps per event 
#   o  information on the number of tracks per event
#   o  information on track length
#   o  information on particle exiting the calorimeter
#
# See  ***LOOKHERE***  to switch on/off some of this observables.
#
#-------------------------------------------------------------------

import os
import sys
import string
import math

#***LOOKHERE***
isShowerPerParticleInformationOn = 1
isNumberOfStepsInformationOn = 1 
isNumberOfTracksInformationOn = 1
isTrackLengthInformationOn = 1
isExitingInformationOn = 1


#===============================================
#================= FUNCTIONS =================== 
#===============================================

def funExtract( theFile ) :
    # Given the file in input, it returns all the observables values,
    # in the form of pairs: value error .
    
    if ( theFile ) :
        statusEndRun = 0
        statusLfractions = 0
        statusRfractions = 0
        statusShowerPerParticle = 0
        statusShowerElectron = 0
        statusShowerProton = 0
        statusShowerPion = 0
        statusShowerNucleus = 0
        statusShowerMuon = 0
        statusShowerKaon = 0
        statusSteps = 0
        statusTracks = 0
        statusTrackLengths = 0
        statusStepLengths = 0
        statusExiting = 0
        for line in theFile :
            #print line
            if ( line.find( "Run Summary" ) > -1 ) :
                statusEndRun = 1
            if ( statusEndRun ) :
                if ( line.find( "longitudinal fraction" ) > -1 ) :
                    statusLfractions = 1
                elif( line.find( "transverse fraction" ) > -1 ) :
                    statusRfractions = 1
                    statusLfractions = 0
                elif ( line.find( "Contributions of the main particle types" ) > -1 ) :
                    statusShowerPerParticle = 1
                    statusRfractions = 0
                elif ( line.find( "STEPS" ) > -1 ) :
                    statusSteps = 1
                    statusShowerPerParticle = 0
                elif ( line.find( "TRACKS" ) > -1 ) :
                    statusTracks = 1
                    statusSteps = 0
                elif ( line.find( "Average track LENGTH" ) > -1 ) :
                    statusTrackLengths = 1
                    statusTracks = 0
                elif ( line.find( "steps and step LENGTH" ) > -1 ) :
                    statusStepLengths = 1
                    statusTrackLengths = 0
                elif ( line.find( "Average exiting" ) > -1 ) :
                    statusExiting = 1
                    statusStepLengths = 0
                    
                if ( not statusShowerPerParticle ) :
                    if ( line.find( "User=" ) > -1 ) :
                        print "  cpuTime = ", \
                              float( line.split( "s " )[0].split( "=" )[1] )
                    elif ( line.find( "active layers =" ) > -1 ) :
                        print "  visEnergy = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "whole calorimeter =" ) > -1 ) :
                        print "  totEnergy = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "energy resolution =" ) > -1 ) : 
                        resolutionValue = float( line.split( "=" )[1].split( "+/-" )[0] )
                        resolutionValue *= 100.0
                        print "  resolution = %.2f" % resolutionValue

                    if ( statusLfractions  and  not statusRfractions ) :
                        if ( line.find( "1st" ) > -1 ) :
                            print "  fL1 = %.1f" % \
                                  float( line.split( "=" )[1].split( "+/-" )[0] )
                        elif ( line.find( "2nd" ) > -1 ) :
                            print "  fL2 = %.1f" % \
                                  float( line.split( "=" )[1].split( "+/-" )[0] )
                        elif ( line.find( "3rd" ) > -1 ) :
                            print "  fL3 = %.1f" % \
                                  float( line.split( "=" )[1].split( "+/-" )[0] )
                        elif ( line.find( "4th" ) > -1 ) :
                            print "  fL4 = %.1f" % \
                                  float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( not statusLfractions  and statusRfractions ) :
                        if ( line.find( "1st" ) > -1 ) :
                            print "  fR1 = %.1f" % \
                                  float( line.split( "=" )[1].split( "+/-" )[0] )
                        elif ( line.find( "2nd" ) > -1 ) :
                            print "  fR2 = %.1f" % \
                                  float( line.split( "=" )[1].split( "+/-" )[0] )
                        elif ( line.find( "3rd" ) > -1 ) :
                            print "  fR3 = %.1f" % \
                                  float( line.split( "=" )[1].split( "+/-" )[0] )
                            
                else :
                    if ( line.find( "Particle type:" ) > -1 ) :
                        statusRfractions = 0
                        if ( line.find( "electron" ) > -1 ) :
                            statusShowerElectron = 1
                        elif ( line.find( "muon" ) > -1 ) :
                            statusShowerMuon = 1
                            statusShowerElectron = 0
                        elif ( line.find( "pion" ) > -1 ) :
                            statusShowerPion = 1
                            statusShowerMuon = 0
                        elif ( line.find( "kaon" ) > -1 ) :
                            statusShowerKaon = 1
                            statusShowerPion = 0
                        elif ( line.find( "proton" ) > -1 ) :
                            statusShowerProton = 1
                            statusShowerKaon = 0
                        elif ( line.find( "pdg0" ) > -1  or     # Before G4 8.2 nuclei 
                               line.find( "nuclei" ) > -1       # have PDG code = 0 .
                               ) :
                            statusShowerNucleus = 1
                            statusShowerProton = 0

                    if ( statusShowerElectron  and  isShowerPerParticleInformationOn ) :
                        if ( line.find( "<E_vis> =" ) > -1 ) :
                            print "  electronEvis = %.1f" % \
                                  float( line.split( "(" )[1].split( "+/-" )[0] )
                        elif ( line.find( "<E_tot> =" ) > -1 ) :
                            print "  electronEtot = %.1f" % \
                                  float( line.split( "(" )[1].split( "+/-" )[0] )

                        if ( statusLfractions  and  not statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                print "  electronfL1 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "2nd" ) > -1 ) :
                                print "  electronfL2 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "3rd" ) > -1 ) :
                                print "  electronfL3 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "4th" ) > -1 ) :
                                print "  electronfL4 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                        elif ( not statusLfractions  and statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                print "  electronfR1 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "2nd" ) > -1 ) :
                                print "  electronfR2 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "3rd" ) > -1 ) :
                                print "  electronfR3 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                                
                    elif ( statusShowerProton  and  isShowerPerParticleInformationOn ) :
                        if ( line.find( "<E_vis> =" ) > -1 ) :
                            print "  protonEvis = %.1f" % \
                                  float( line.split( "(" )[1].split( "+/-" )[0] )
                        elif ( line.find( "<E_tot> =" ) > -1 ) :
                            print "  protonEtot = %.1f" % \
                                  float( line.split( "(" )[1].split( "+/-" )[0] )

                        if ( statusLfractions  and  not statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                print "  protonfL1 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "2nd" ) > -1 ) :
                                print "  protonfL2 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "3rd" ) > -1 ) :
                                print "  protonfL3 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "4th" ) > -1 ) :
                                print "  protonfL4 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                        elif ( not statusLfractions  and statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                print "  protonfR1 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "2nd" ) > -1 ) :
                                print "  protonfR2 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "3rd" ) > -1 ) :
                                print "  protonfR3 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                                
                    elif ( statusShowerPion  and  isShowerPerParticleInformationOn ) :
                        if ( line.find( "<E_vis> =" ) > -1 ) :
                            print "  pionEvis = %.1f" % \
                                  float( line.split( "(" )[1].split( "+/-" )[0] )
                        elif ( line.find( "<E_tot> =" ) > -1 ) :
                            print "  pionEtot = %.1f" % \
                                  float( line.split( "(" )[1].split( "+/-" )[0] )

                        if ( statusLfractions  and  not statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                print "  pionfL1 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "2nd" ) > -1 ) :
                                print "  pionfL2 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "3rd" ) > -1 ) :
                                print "  pionfL3 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "4th" ) > -1 ) :
                                print "  pionfL4 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                        elif ( not statusLfractions  and statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                print "  pionfR1 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "2nd" ) > -1 ) :
                                print "  pionfR2 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "3rd" ) > -1 ) :
                                print "  pionfR3 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                                
                    elif ( statusShowerNucleus  and  isShowerPerParticleInformationOn ) :
                        if ( line.find( "<E_vis> =" ) > -1 ) :
                            print "  nucleusEvis = %.1f" % \
                                  float( line.split( "(" )[1].split( "+/-" )[0] )
                        elif ( line.find( "<E_tot> =" ) > -1 ) :
                            print "  nucleusEtot = %.1f" % \
                                  float( line.split( "(" )[1].split( "+/-" )[0] )

                        if ( statusLfractions  and  not statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                print "  nucleusfL1 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "2nd" ) > -1 ) :
                                print "  nucleusfL2 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "3rd" ) > -1 ) :
                                print "  nucleusfL3 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "4th" ) > -1 ) :
                                print "  nucleusfL4 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                        elif ( not statusLfractions  and statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                print "  nucleusfR1 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "2nd" ) > -1 ) :
                                print "  nucleusfR2 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )
                            elif ( line.find( "3rd" ) > -1 ) :
                                print "  nucleusfR3 = %.1f" % \
                                      float( line.split( "=" )[1].split( "+/-" )[0] )

                    elif ( statusShowerMuon  and  isShowerPerParticleInformationOn ) :
                        if ( line.find( "<E_vis> =" ) > -1 ) :
                            print "  muonEvis = %.1f" % \
                                  float( line.split( "(" )[1].split( "+/-" )[0] )
                        elif ( line.find( "<E_tot> =" ) > -1 ) :
                            print "  muonEtot = %.1f" % \
                                  float( line.split( "(" )[1].split( "+/-" )[0] )

                    elif ( statusShowerKaon  and  isShowerPerParticleInformationOn ) :
                        if ( line.find( "<E_vis> =" ) > -1 ) :
                            print "  kaonEvis = %.1f" % \
                                  float( line.split( "(" )[1].split( "+/-" )[0] )
                        elif ( line.find( "<E_tot> =" ) > -1 ) :
                            print "  kaonEtot = %.1f" % \
                                  float( line.split( "(" )[1].split( "+/-" )[0] )

                if ( statusSteps  and  isNumberOfStepsInformationOn ) :
                    if ( line.find( "+/-" ) > -1 ) :
                        print line.split( "=" )[0], "steps = %.0f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )

                if ( statusTracks  and  isNumberOfTracksInformationOn ) :
                    if ( line.find( "+/-" ) > -1 ) :
                        print line.split( "=" )[0], "tracks = %.0f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )

                if ( statusTrackLengths  and  isTrackLengthInformationOn ) :
                    if ( line.find( "electron/positron" ) > -1 ) :
                        print "  electronLength = %.2f" % \
                              float( line.split( ":" )[1].split( "+/-" )[0] )
                    elif ( line.find( "pion-/pion+" ) > -1 ) :
                        print "  pionLength = %.1f" % \
                              float( line.split( ":" )[1].split( "+/-" )[0] )
                    elif ( line.find( "proton" ) > -1 ) :
                        print "  protonLength = %.1f" % \
                              float( line.split( ":" )[1].split( "+/-" )[0] )
                    elif ( line.find( "gamma" ) > -1 ) :
                        print "  gammaLength = %.1f" % \
                              float( line.split( ":" )[1].split( "+/-" )[0] )
                    elif ( line.find( "neutron" ) > -1 ) :
                        print "  neutronLength = %.1f" % \
                              float( line.split( ":" )[1].split( "+/-" )[0] )
                        
                if ( statusExiting  and  isExitingInformationOn ) :
                    if ( line.find( "exiting Kinetic Energy" ) > -1 ) :
                        print "  exitKin = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "fraction due to Gammas" ) > -1 ) :
                        print "  exitFracGammas = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "fraction due to Neutrons" ) > -1 ) :
                        print "  exitFracNeutrons = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "fraction due to Neutrinos" ) > -1 ) :
                        print "  exitFracNeutrinos = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "fraction due to Muons" ) > -1 ) :
                        print "  exitFracMuons = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "fraction due to Electrons" ) > -1 ) :
                        print "  exitFracElectrons = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "fraction due to Others" ) > -1 ) :
                        print "  exitFracOthers = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "exiting particles" ) > -1 ) :
                        print "  exitNum = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "exiting Gammas" ) > -1 ) :
                        print "  exitNumGammas = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "exiting Neutrons" ) > -1 ) :
                        print "  exitNumNeutrons = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "exiting Neutrinos" ) > -1 ) :
                        print "  exitNumNeutrinos = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "exiting Muons" ) > -1 ) :
                        print "  exitNumMuons = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "exiting Electrons" ) > -1 ) :
                        print "  exitNumElectrons = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )
                    elif ( line.find( "exiting Others" ) > -1 ) :
                        print "  exitNumOthers = %.1f" % \
                              float( line.split( "=" )[1].split( "+/-" )[0] )

    return


#===============================================
#==================== MAIN ===================== 
#===============================================

# Check whether the number of arguments are correct or not.
if ( len( sys.argv ) != 2 ) :
    print " Usage:  printInfoLogfile.py file"
else :
    theFile = open( sys.argv[1], "r" )
    funExtract( theFile )
    theFile.close()

#-------------------------------------------------------------------

