#!/usr/bin/python

#-------------------------------------------------------------------
# Last update: 15-Jul-2008
#
# This script has 1 input argument, and should be run as:
#
#         $  python printInfoLogfile.py file
#
# where the input  file  should be one of the log-files
# produced by running the StatAccepTest simulations.
# This script prints out (in the screen) the following
# observables:
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
# as extracted from the log-files.
# In the case that the  file  does not exist, or it exists
# but is not a log-file obtained by completing a running
# of StatAccepTest, then "0" is print for each of the above
# observables.
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
    
    statusEndRun = 0
    if ( theFile ) :
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

    if ( statusEndRun == 0 ) :
        print "  cpuTime =  0  "
        print "  visEnergy =  0  "
        print "  totEnergy =  0  "
        print "  resolution =  0  "
        print "  fL1 =  0  "
        print "  fL2 =  0  "
        print "  fL3 =  0  "
        print "  fL4 =  0  "
        print "  fR1 =  0  "
        print "  fR2 =  0  "
        print "  fR3 =  0  "
        if ( isShowerPerParticleInformationOn ) :
            print "  electronEvis =  0  "
            print "  electronEtot =  0  "
            print "  electronfL1 =  0  "
            print "  electronfL2 =  0  "
            print "  electronfL3 =  0  "
            print "  electronfL4 =  0  "
            print "  electronfR1 =  0  "
            print "  electronfR2 =  0  "
            print "  electronfR3 =  0  "
            print "  protonEvis =  0  "
            print "  protonEtot =  0  "
            print "  protonfL1 =  0  "
            print "  protonfL2 =  0  "
            print "  protonfL3 =  0  "
            print "  protonfL4 =  0  "
            print "  protonfR1 =  0  "
            print "  protonfR2 =  0  "
            print "  protonfR3 =  0  "
            print "  pionEvis =  0  "
            print "  pionEtot =  0  "
            print "  pionfL1 =  0  "
            print "  pionfL2 =  0  "
            print "  pionfL3 =  0  "
            print "  pionfL4 =  0  "
            print "  pionfR1 =  0  "
            print "  pionfR2 =  0  "
            print "  pionfR3 =  0  "
            print "  nucleusEvis =  0  "
            print "  nucleusEtot =  0  "
            print "  nucleusfL1 =  0  "
            print "  nucleusfL2 =  0  "
            print "  nucleusfL3 =  0  "
            print "  nucleusfL4 =  0  "
            print "  nucleusfR1 =  0  "
            print "  nucleusfR2 =  0  "
            print "  nucleusfR3 =  0  "
            print "  muonEvis =  0  "
            print "  muonEtot =  0  "
            print "  kaonEvis =  0  "
            print "  kaonEtot =  0  "
        if ( isNumberOfStepsInformationOn ) :
            print "  # total steps =  0  " 
            print "  # positives steps =  0  " 
            print "  # neutrals steps =  0  " 
            print "  # negatives steps =  0  " 
            print "  # nuclei steps =  0  " 
            print "  # particles with Unrecognized PDG code steps =  0  " 
            print "  # electromagnetic (e+ , e- , gammas) steps =  0  " 
            print "  # electroweak (mu+, mu-, tau+, tau-, neutrinos) steps =  0  " 
            print "  # hadrons steps =  0  " 
            print "  # mesons steps =  0  " 
            print "  # baryons steps =  0  " 
            print "  # light mesons (u/ubar/d/dbar) steps =  0  " 
            print "  # light baryons (u/ubar/d/dbar) steps =  0  " 
            print "  # strange (s/sbar) mesons steps =  0  " 
            print "  # strange (s/sbar) baryons steps =  0  " 
            print "  # heavy (c/cbar or b/bbar) mesons steps =  0  " 
            print "  # heavy (c/cbar or b/bbar) baryons steps =  0  " 
            print "  # electrons steps =  0  " 
            print "  # gammas steps =  0  " 
            print "  # positrons steps =  0  " 
            print "  # mu- steps =  0  " 
            print "  # mu+ steps =  0  " 
            print "  # tau- steps =  0  "
            print "  # tau+ steps =  0  " 
            print "  # neutrinos steps =  0  " 
            print "  # pi+ steps =  0  " 
            print "  # pi0 steps =  0  " 
            print "  # pi- steps =  0  " 
            print "  # K+ steps =  0  " 
            print "  # K-neutral (K0/K0bar or K0_S/K0_L) steps =  0  " 
            print "  # K- steps =  0  " 
            print "  # protons steps =  0  " 
            print "  # anti-protons steps =  0  " 
            print "  # neutrons steps =  0  " 
            print "  # anti-neutrons steps =  0  " 
        if ( isNumberOfTracksInformationOn ) :
            print "  # total tracks =  0  " 
            print "  # positives tracks =  0  " 
            print "  # neutrals tracks =  0  " 
            print "  # negatives tracks =  0  " 
            print "  # nuclei tracks =  0  " 
            print "  # particles with Unrecognized PDG code tracks =  0  " 
            print "  # electromagnetic (e+ , e- , gammas) tracks =  0  " 
            print "  # electroweak (mu+, mu-, tau+, tau-, neutrinos) tracks =  0  " 
            print "  # hadrons tracks =  0  " 
            print "  # mesons tracks =  0  " 
            print "  # baryons tracks =  0  " 
            print "  # light mesons (u/ubar/d/dbar) tracks =  0  " 
            print "  # light baryons (u/ubar/d/dbar) tracks =  0  " 
            print "  # strange (s/sbar) mesons tracks =  0  " 
            print "  # strange (s/sbar) baryons tracks =  0  " 
            print "  # heavy (c/cbar or b/bbar) mesons tracks =  0  " 
            print "  # heavy (c/cbar or b/bbar) baryons tracks =  0  " 
            print "  # electrons tracks =  0  " 
            print "  # gammas tracks =  0  " 
            print "  # positrons tracks =  0  " 
            print "  # mu- tracks =  0  " 
            print "  # mu+ tracks =  0  " 
            print "  # tau- tracks =  0  "
            print "  # tau+ tracks =  0  " 
            print "  # neutrinos tracks =  0  " 
            print "  # pi+ tracks =  0  " 
            print "  # pi0 tracks =  0  " 
            print "  # pi- tracks =  0  " 
            print "  # K+ tracks =  0  " 
            print "  # K-neutral (K0/K0bar or K0_S/K0_L) tracks =  0  " 
            print "  # K- tracks =  0  " 
            print "  # protons tracks =  0  " 
            print "  # anti-protons tracks =  0  " 
            print "  # neutrons tracks =  0  " 
            print "  # anti-neutrons tracks =  0  "   
        if ( isTrackLengthInformationOn ) :
            print "  electronLength =  0  "
            print "  pionLength =  0  "
            print "  protonLength =  0  "
            print "  gammaLength =  0  "
            print "  neutronLength =  0  "
        if ( isExitingInformationOn ) :
            print "  exitKin =  0  "
            print "  exitFracGammas =  0  "
            print "  exitFracNeutrons =  0  "
            print "  exitFracNeutrinos =  0  "
            print "  exitFracMuons =  0  "
            print "  exitFracElectrons =  0  "
            print "  exitFracOthers =  0  "
            print "  exitNum =  0  "
            print "  exitNumGammas =  0  "
            print "  exitNumNeutrons =  0  "
            print "  exitNumNeutrinos =  0  "
            print "  exitNumMuons =  0  "
            print "  exitNumElectrons =  0  "
            print "  exitNumOthers =  0  "
            
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

