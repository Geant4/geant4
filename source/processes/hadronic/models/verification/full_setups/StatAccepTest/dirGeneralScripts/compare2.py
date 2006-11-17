#!/usr/bin/python

#-------------------------------------------------------------------
# Last update: 17-Nov-2006
#
# This script has 2 input arguments, and should be run as:
#
#         $  python compare2 fileA fileB
#
# where the two input files should be the log file of running
# the StatAccepTest simulations.
# This script compares the following observables between the
# two files, and print them in case that they differ by more
# than  5 sigma  and/or  by more than  5%
# ( NB) See  ***LOOKHERE***  for the choice of the thresholds
#       and whether the condition is based on "and" or "or".)
#
#   o  CPU times (without error)
#   o  total visible energy
#   o  total energy in the whole calorimeter
#   o  visible energy in each layer (in MeV)
#   o  fractions in 4 longitudinal sections
#   o  visible energy in each ring (in MeV)
#   o  fractions in 3 radial sections 
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
thresholdSigmaDifference = 5.0
thresholdRelativeDifference = 0.05
isConditionAND = 0    # 1 = AND  ; 0 = OR

isLayerInformationOn = 0
isRingInformationOn = 0
isShowerPerParticleInformationOn = 1
isNumberOfStepsInformationOn = 0
isNumberOfTracksInformationOn = 1
isTrackLengthInformationOn = 1
isExitingInformationOn = 1


#===============================================
#================= FUNCTIONS =================== 
#===============================================

def funExtract( theFile ) :
    # Given the file in input, it returns all the observables values,
    # in the form of pairs: (value, error).
    # In the case of collections of related values, as the visible
    # energy in each layer, the visible energy in each ring, and the
    # information on the average number of steps and tracks per event,
    # the pairs are grouped in lists.
    
    cpuTime = 0.0
    visEnergy = (0.0, 0.0)
    totEnergy = (0.0, 0.0)
    listL = []
    listR = []
    listSteps = []
    listTracks = []
    fL1 = (0.0, 0.0)
    fL2 = (0.0, 0.0)
    fL3 = (0.0, 0.0)
    fL4 = (0.0, 0.0)
    fR1 = (0.0, 0.0)
    fR2 = (0.0, 0.0)
    fR3 = (0.0, 0.0)
    electronEvis = (0.0, 0.0)
    electronEtot = (0.0, 0.0)
    electronfL1 = (0.0, 0.0)
    electronfL2 = (0.0, 0.0)
    electronfL3 = (0.0, 0.0)
    electronfL4 = (0.0, 0.0)
    electronfR1 = (0.0, 0.0)
    electronfR2 = (0.0, 0.0)
    electronfR3 = (0.0, 0.0)
    protonEvis = (0.0, 0.0)
    protonEtot = (0.0, 0.0)
    protonfL1 = (0.0, 0.0)
    protonfL2 = (0.0, 0.0)
    protonfL3 = (0.0, 0.0)
    protonfL4 = (0.0, 0.0)
    protonfR1 = (0.0, 0.0)
    protonfR2 = (0.0, 0.0)
    protonfR3 = (0.0, 0.0)
    pionEvis = (0.0, 0.0)
    pionEtot = (0.0, 0.0)
    pionfL1 = (0.0, 0.0)
    pionfL2 = (0.0, 0.0)
    pionfL3 = (0.0, 0.0)
    pionfL4 = (0.0, 0.0)
    pionfR1 = (0.0, 0.0)
    pionfR2 = (0.0, 0.0)
    pionfR3 = (0.0, 0.0)
    nucleusEvis = (0.0, 0.0)
    nucleusEtot = (0.0, 0.0)
    nucleusfL1 = (0.0, 0.0)
    nucleusfL2 = (0.0, 0.0)
    nucleusfL3 = (0.0, 0.0)
    nucleusfL4 = (0.0, 0.0)
    nucleusfR1 = (0.0, 0.0)
    nucleusfR2 = (0.0, 0.0)
    nucleusfR3 = (0.0, 0.0)
    muonEvis = (0.0, 0.0)
    muonEtot = (0.0, 0.0)
    kaonEvis = (0,0, 0.0)
    kaonEtot = (0,0, 0.0)
    electronLength = (0.0, 0.0)
    pionLength = (0.0, 0.0)
    protonLength = (0.0, 0.0)
    gammaLength = (0.0, 0.0)
    neutronLength = (0.0, 0.0)
    exitKin = (0.0, 0.0)
    exitFracGammas = (0.0, 0.0)
    exitFracNeutrons = (0.0, 0.0)
    exitFracNeutrinos = (0.0, 0.0)
    exitFracMuons = (0.0, 0.0)
    exitFracElectrons = (0.0, 0.0)
    exitFracOthers = (0.0, 0.0)
    exitNum = (0.0, 0.0)
    exitNumGammas = (0.0, 0.0)
    exitNumNeutrons = (0.0, 0.0)
    exitNumNeutrinos = (0.0, 0.0)
    exitNumMuons = (0.0, 0.0)
    exitNumElectrons = (0.0, 0.0)
    exitNumOthers = (0.0, 0.0)

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
                        cpuTime = float( line.split( "s " )[0].split( "=" )[1] )
                        #print "***DEBUG***  cpuTime = ", cpuTime
                    elif ( line.find( "active layers =" ) > -1 ) :
                        visEnergy = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                      float( line.split( "=" )[1].split( "+/-" )[1] ) )
                        #print "***DEBUG***  visEnergy = ", visEnergy
                    elif ( line.find( "whole calorimeter =" ) > -1 ) :
                        totEnergy = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                      float( line.split( "=" )[1].split( "+/-" )[1] ) )
                        #print "***DEBUG***  totEnergy = ", totEnergy
                    elif ( line.find( "layer =" ) > -1 ) :
                        listL.append( ( float( line.split( "=" )[2].split( "+/-" )[0] ),
                                        float( line.split( "=" )[2].split( "+/-" )[1] ) ) )
                        #print "***DEBUG***  listL = ", listL
                    elif ( line.find( "iBinR =" ) > -1 ) :
                        listR.append( ( float( line.split( "=" )[2].split( "+/-" )[0] ),
                                        float( line.split( "=" )[2].split( "+/-" )[1] ) ) )
                        #print "***DEBUG***  listR = ", listR

                    if ( statusLfractions  and  not statusRfractions ) :
                        if ( line.find( "1st" ) > -1 ) :
                            fL1 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                    float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                            #print "***DEBUG***  fL1 = ", fL1
                        elif ( line.find( "2nd" ) > -1 ) :
                            fL2 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                    float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                            #print "***DEBUG***  fL2 = ", fL2
                        elif ( line.find( "3rd" ) > -1 ) :
                            fL3 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                    float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                            #print "***DEBUG***  fL3 = ", fL3
                        elif ( line.find( "4th" ) > -1 ) :
                            fL4 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                    float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                            #print "***DEBUG***  fL4 = ", fL4
                    elif ( not statusLfractions  and statusRfractions ) :
                        if ( line.find( "1st" ) > -1 ) :
                            fR1 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                    float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                            #print "***DEBUG***  fR1 = ", fR1
                        elif ( line.find( "2nd" ) > -1 ) :
                            fR2 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                    float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                            #print "***DEBUG***  fR2 = ", fR2
                        elif ( line.find( "3rd" ) > -1 ) :
                            fR3 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                    float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                            #print "***DEBUG***  fR3 = ", fR3
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

                    if ( statusShowerElectron ) :
                        if ( line.find( "<E_vis> =" ) > -1 ) :
                            electronEvis = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                             float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                            #print "***DEBUG***  electronEvis = ", electronEvis
                        elif ( line.find( "<E_tot> =" ) > -1 ) :
                            electronEtot = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                             float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                            #print "***DEBUG***  electronEtot = ", electronEtot

                        if ( statusLfractions  and  not statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                electronfL1 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                                float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  electronfL1 = ", electronfL1
                            elif ( line.find( "2nd" ) > -1 ) :
                                electronfL2 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                                float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  electronfL2 = ", electronfL2
                            elif ( line.find( "3rd" ) > -1 ) :
                                electronfL3 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                                float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  electronfL3 = ", electronfL3
                            elif ( line.find( "4th" ) > -1 ) :
                                electronfL4 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                                float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  electronfL4 = ", electronfL4
                        elif ( not statusLfractions  and statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                electronfR1 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                                float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  electronfR1 = ", electronfR1
                            elif ( line.find( "2nd" ) > -1 ) :
                                electronfR2 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                                float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  electronfR2 = ", electronfR2
                            elif ( line.find( "3rd" ) > -1 ) :
                                electronfR3 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                                float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  electronfR3 = ", electronfR3
                                
                    elif ( statusShowerProton ) :
                        if ( line.find( "<E_vis> =" ) > -1 ) :
                            protonEvis = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                           float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                            #print "***DEBUG***  protonEvis = ", protonEvis
                        elif ( line.find( "<E_tot> =" ) > -1 ) :
                            protonEtot = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                           float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                            #print "***DEBUG***  protonEtot = ", protonEtot

                        if ( statusLfractions  and  not statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                protonfL1 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                              float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  protonfL1 = ", protonfL1
                            elif ( line.find( "2nd" ) > -1 ) :
                                protonfL2 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                              float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  protonfL2 = ", protonfL2
                            elif ( line.find( "3rd" ) > -1 ) :
                                protonfL3 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                              float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  protonfL3 = ", protonfL3
                            elif ( line.find( "4th" ) > -1 ) :
                                protonfL4 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                              float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  protonfL4 = ", protonfL4
                        elif ( not statusLfractions  and statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                protonfR1 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                              float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  protonfR1 = ", protonfR1
                            elif ( line.find( "2nd" ) > -1 ) :
                                protonfR2 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                              float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  protonfR2 = ", protonfR2
                            elif ( line.find( "3rd" ) > -1 ) :
                                protonfR3 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                              float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  protonfR3 = ", protonfR3
                                
                    elif ( statusShowerPion ) :
                        if ( line.find( "<E_vis> =" ) > -1 ) :
                            pionEvis = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                         float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                            #print "***DEBUG***  pionEvis = ", pionEvis
                        elif ( line.find( "<E_tot> =" ) > -1 ) :
                            pionEtot = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                         float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                            #print "***DEBUG***  pionEtot = ", pionEtot

                        if ( statusLfractions  and  not statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                pionfL1 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  pionfL1 = ", pionfL1
                            elif ( line.find( "2nd" ) > -1 ) :
                                pionfL2 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  pionfL2 = ", pionfL2
                            elif ( line.find( "3rd" ) > -1 ) :
                                pionfL3 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  pionfL3 = ", pionfL3
                            elif ( line.find( "4th" ) > -1 ) :
                                pionfL4 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  pionfL4 = ", pionfL4
                        elif ( not statusLfractions  and statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                pionfR1 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  pionfR1 = ", pionfR1
                            elif ( line.find( "2nd" ) > -1 ) :
                                pionfR2 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  pionfR2 = ", pionfR2
                            elif ( line.find( "3rd" ) > -1 ) :
                                pionfR3 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  pionfR3 = ", pionfR3
                                
                    elif ( statusShowerNucleus ) :
                        if ( line.find( "<E_vis> =" ) > -1 ) :
                            nucleusEvis = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                         float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                            #print "***DEBUG***  nucleusEvis = ", nucleusEvis
                        elif ( line.find( "<E_tot> =" ) > -1 ) :
                            nucleusEtot = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                         float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                            #print "***DEBUG***  nucleusEtot = ", nucleusEtot

                        if ( statusLfractions  and  not statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                nucleusfL1 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  nucleusfL1 = ", nucleusfL1
                            elif ( line.find( "2nd" ) > -1 ) :
                                nucleusfL2 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  nucleusfL2 = ", nucleusfL2
                            elif ( line.find( "3rd" ) > -1 ) :
                                nucleusfL3 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  nucleusfL3 = ", nucleusfL3
                            elif ( line.find( "4th" ) > -1 ) :
                                nucleusfL4 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  nucleusfL4 = ", nucleusfL4
                        elif ( not statusLfractions  and statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                nucleusfR1 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  nucleusfR1 = ", nucleusfR1
                            elif ( line.find( "2nd" ) > -1 ) :
                                nucleusfR2 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  nucleusfR2 = ", nucleusfR2
                            elif ( line.find( "3rd" ) > -1 ) :
                                nucleusfR3 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  nucleusfR3 = ", nucleusfR3

                    elif ( statusShowerMuon ) :
                        if ( line.find( "<E_vis> =" ) > -1 ) :
                            muonEvis = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                         float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                            #print "***DEBUG***  muonEvis = ", muonEvis
                        elif ( line.find( "<E_tot> =" ) > -1 ) :
                            muonEtot = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                         float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                            #print "***DEBUG***  muonEtot = ", muonEtot

                    elif ( statusShowerKaon ) :
                        if ( line.find( "<E_vis> =" ) > -1 ) :
                            kaonEvis = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                         float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                            #print "***DEBUG***  kaonEvis = ", kaonEvis
                        elif ( line.find( "<E_tot> =" ) > -1 ) :
                            kaonEtot = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                         float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                            #print "***DEBUG***  kaonEtot = ", kaonEtot

                if ( statusSteps ) :
                    if ( line.find( "+/-" ) > -1 ) :
                        listSteps.append( ( float( line.split( "=" )[1].split( "+/-" )[0] ) , float( line.split( "=" )[1].split( "+/-" )[1] ) ) )

                if ( statusTracks ) :
                    if ( line.find( "+/-" ) > -1 ) :
                        listTracks.append( ( float( line.split( "=" )[1].split( "+/-" )[0] ) , float( line.split( "=" )[1].split( "+/-" )[1] ) ) )

                if ( statusTrackLengths ) :
                    if ( line.find( "electron/positron" ) > -1 ) :
                        electronLength = ( float( line.split( ":" )[1].split( "+/-" )[0] ),
                                           float( line.split( ":" )[1].split( "+/-" )[1] ) )
                        #print "***DEBUG***  electronLength = ", electronLength
                    elif ( line.find( "pion-/pion+" ) > -1 ) :
                        pionLength = ( float( line.split( ":" )[1].split( "+/-" )[0] ),
                                       float( line.split( ":" )[1].split( "+/-" )[1] ) )
                        #print "***DEBUG***  pionLength = ", pionLength
                    elif ( line.find( "proton" ) > -1 ) :
                        protonLength = ( float( line.split( ":" )[1].split( "+/-" )[0] ),
                                         float( line.split( ":" )[1].split( "+/-" )[1] ) )
                        #print "***DEBUG***  protonLength = ", protonLength
                    elif ( line.find( "gamma" ) > -1 ) :
                        gammaLength = ( float( line.split( ":" )[1].split( "+/-" )[0] ),
                                        float( line.split( ":" )[1].split( "+/-" )[1] ) )
                        #print "***DEBUG***  gammaLength = ", gammaLength
                    elif ( line.find( "neutron" ) > -1 ) :
                        neutronLength = ( float( line.split( ":" )[1].split( "+/-" )[0] ),
                                          float( line.split( ":" )[1].split( "+/-" )[1] ) )
                        #print "***DEBUG***  neutronLength = ", neutronLength

                if ( statusExiting ) :
                    if ( line.find( "exiting Kinetic Energy" ) > -1 ) :
                        exitKin = ( float( line.split( "=" )[1].split( "+/-" )[0] ),
                                    float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                        #print "***DEBUG***  exitKin = ", exitKin
                    elif ( line.find( "fraction due to Gammas" ) > -1 ) :
                        exitFracGammas = ( float( line.split( "=" )[1].split( "+/-" )[0] ),
                                           float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                        #print "***DEBUG***  exitFracGammas = ", exitFracGammas
                    elif ( line.find( "fraction due to Neutrons" ) > -1 ) :
                        exitFracNeutrons = ( float( line.split( "=" )[1].split( "+/-" )[0] ),
                                             float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                        #print "***DEBUG***  exitFracNeutrons = ", exitFracNeutrons
                    elif ( line.find( "fraction due to Neutrinos" ) > -1 ) :
                        exitFracNeutrinos = ( float( line.split( "=" )[1].split( "+/-" )[0] ),
                                              float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                        #print "***DEBUG***  exitFracNeutrinos = ", exitFracNeutrinos
                    elif ( line.find( "fraction due to Muons" ) > -1 ) :
                        exitFracMuons = ( float( line.split( "=" )[1].split( "+/-" )[0] ),
                                          float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                        #print "***DEBUG***  exitFracMuons = ", exitFracMuons
                    elif ( line.find( "fraction due to Electrons" ) > -1 ) :
                        exitFracElectrons = ( float( line.split( "=" )[1].split( "+/-" )[0] ),
                                              float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                        #print "***DEBUG***  exitFracElectrons = ", exitFracElectrons
                    elif ( line.find( "fraction due to Others" ) > -1 ) :
                        exitFracOthers = ( float( line.split( "=" )[1].split( "+/-" )[0] ),
                                           float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                        #print "***DEBUG***  exitFracOthers = ", exitFracOthers
                    elif ( line.find( "exiting particles" ) > -1 ) :
                        exitNum = ( float( line.split( "=" )[1].split( "+/-" )[0] ),
                                    float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                        #print "***DEBUG***  exitNum = ", exitNum
                    elif ( line.find( "exiting Gammas" ) > -1 ) :
                        exitNumGammas = ( float( line.split( "=" )[1].split( "+/-" )[0] ),
                                          float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                        #print "***DEBUG***  exitNumGammas = ", exitNumGammas
                    elif ( line.find( "exiting Neutrons" ) > -1 ) :
                        exitNumNeutrons = ( float( line.split( "=" )[1].split( "+/-" )[0] ),
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                        #print "***DEBUG***  exitNumNeutrons = ", exitNumNeutrons
                    elif ( line.find( "exiting Neutrinos" ) > -1 ) :
                        exitNumNeutrinos = ( float( line.split( "=" )[1].split( "+/-" )[0] ),
                                             float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                        #print "***DEBUG***  exitNumNeutrinos = ", exitNumNeutrinos
                    elif ( line.find( "exiting Muons" ) > -1 ) :
                        exitNumMuons = ( float( line.split( "=" )[1].split( "+/-" )[0] ),
                                         float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                        #print "***DEBUG***  exitNumMuons = ", exitNumMuons
                    elif ( line.find( "exiting Electrons" ) > -1 ) :
                        exitNumElectrons = ( float( line.split( "=" )[1].split( "+/-" )[0] ),
                                             float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                        #print "***DEBUG***  exitNumElectrons = ", exitNumElectrons
                    elif ( line.find( "exiting Others" ) > -1 ) :
                        exitNumOthers = ( float( line.split( "=" )[1].split( "+/-" )[0] ),
                                          float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                        #print "***DEBUG***  exitNumOthers = ", exitNumOthers

    return ( cpuTime, visEnergy, totEnergy,
             listL, fL1, fL2, fL3, fL4,
             listR, fR1, fR2, fR3,
             electronEvis, electronEtot,
             electronfL1, electronfL2, electronfL3, electronfL4, 
             electronfR1, electronfR2, electronfR3,
             protonEvis, protonEtot,
             protonfL1, protonfL2, protonfL3, protonfL4, 
             protonfR1, protonfR2, protonfR3,
             pionEvis, pionEtot,
             pionfL1, pionfL2, pionfL3, pionfL4, 
             pionfR1, pionfR2, pionfR3,
             nucleusEvis, nucleusEtot,
             nucleusfL1, nucleusfL2, nucleusfL3, nucleusfL4, 
             nucleusfR1, nucleusfR2, nucleusfR3,
             muonEvis, muonEtot,
             kaonEvis, kaonEtot,
             listSteps, listTracks,
             electronLength, pionLength, protonLength, gammaLength, neutronLength,
             exitKin, exitFracGammas, exitFracNeutrons, exitFracNeutrinos,
             exitFracMuons, exitFracElectrons, exitFracOthers,
             exitNum, exitNumGammas, exitNumNeutrons, exitNumNeutrinos,
             exitNumMuons, exitNumElectrons, exitNumOthers )


def doStatisticalComparison( valuesA, valuesB ) :
    # Given in input the two sets of observables, this function does
    # the comparison of each of these observable, one by one.
    # It computes both the deviations in terms of the error (i.e. the
    # number of sigma of deviation between the two observables), and
    # in terms of relative deviation (i.e. the ratio of the deviation
    # by the average value between the two observables).
    # Only if both these deviations are above a certain threshold
    # (typically 5 sigma and 5% respectively), the two observables
    # are printed out.
    # In the case (as for the CPU time) the observables have no
    # error, then only the relative deviation is considered.

    ( cpuTimeA, visEnergyA, totEnergyA,
      listLA, fL1A, fL2A, fL3A, fL4A,
      listRA, fR1A, fR2A, fR3A,
      electronEvisA, electronEtotA,
      electronfL1A, electronfL2A, electronfL3A, electronfL4A, 
      electronfR1A, electronfR2A, electronfR3A,
      protonEvisA, protonEtotA,
      protonfL1A, protonfL2A, protonfL3A, protonfL4A, 
      protonfR1A, protonfR2A, protonfR3A,
      pionEvisA, pionEtotA,
      pionfL1A, pionfL2A, pionfL3A, pionfL4A, 
      pionfR1A, pionfR2A, pionfR3A,
      nucleusEvisA, nucleusEtotA,
      nucleusfL1A, nucleusfL2A, nucleusfL3A, nucleusfL4A, 
      nucleusfR1A, nucleusfR2A, nucleusfR3A,
      muonEvisA, muonEtotA,
      kaonEvisA, kaonEtotA,
      listStepsA, listTracksA,
      electronLengthA, pionLengthA, protonLengthA, gammaLengthA, neutronLengthA,
      exitKinA, exitFracGammasA, exitFracNeutronsA, exitFracNeutrinosA,
      exitFracMuonsA, exitFracElectronsA, exitFracOthersA,
      exitNumA, exitNumGammasA, exitNumNeutronsA, exitNumNeutrinosA,
      exitNumMuonsA, exitNumElectronsA, exitNumOthersA ) = valuesA
    
    ( cpuTimeB, visEnergyB, totEnergyB,
      listLB, fL1B, fL2B, fL3B, fL4B,
      listRB, fR1B, fR2B, fR3B,
      electronEvisB, electronEtotB,
      electronfL1B, electronfL2B, electronfL3B, electronfL4B, 
      electronfR1B, electronfR2B, electronfR3B,
      protonEvisB, protonEtotB,
      protonfL1B, protonfL2B, protonfL3B, protonfL4B, 
      protonfR1B, protonfR2B, protonfR3B,
      pionEvisB, pionEtotB,
      pionfL1B, pionfL2B, pionfL3B, pionfL4B, 
      pionfR1B, pionfR2B, pionfR3B,
      nucleusEvisB, nucleusEtotB,
      nucleusfL1B, nucleusfL2B, nucleusfL3B, nucleusfL4B, 
      nucleusfR1B, nucleusfR2B, nucleusfR3B,
      muonEvisB, muonEtotB,
      kaonEvisB, kaonEtotB,
      listStepsB, listTracksB,
      electronLengthB, pionLengthB, protonLengthB, gammaLengthB, neutronLengthB,
      exitKinB, exitFracGammasB, exitFracNeutronsB, exitFracNeutrinosB,
      exitFracMuonsB, exitFracElectronsB, exitFracOthersB,
      exitNumB, exitNumGammasB, exitNumNeutronsB, exitNumNeutrinosB,
      exitNumMuonsB, exitNumElectronsB, exitNumOthersB ) = valuesB

    listResults = []  # A list of results, where the result is a tuple:
                      # (nameObservable, valueA, valueB, resultComparison)
                      # where  resultComparison  is, in turn, a tuple:
                      #   ( sigmaDiff, relativeDiff )

    for i in range(73) :
        if ( i == 0 ) :
            name = " CPU times [sec]"
            resultComparison = compareTwo( cpuTimeA, 0.0, cpuTimeB, 0.0 )
            listResults.append( ( name, cpuTimeA, cpuTimeB, resultComparison ) )
        elif ( i == 1 ) :
            name = " Evis [MeV]"
            resultComparison = compareTwo( visEnergyA[0], visEnergyA[1],
                                           visEnergyB[0], visEnergyB[1] )
            listResults.append( ( name, visEnergyA[0], visEnergyB[0], resultComparison ) )
        elif ( i == 2 ) :
            name = " Etot [MeV]"
            resultComparison = compareTwo( totEnergyA[0], totEnergyA[1],
                                           totEnergyB[0], totEnergyB[1] )
            listResults.append( ( name, totEnergyA[0], totEnergyB[0],
                                  resultComparison ) )
        elif ( i == 3 ) :
            if ( len( listLA ) != len( listLB ) ) :
                print " ***DIFFERENT NUMBER OF LAYERS *** : ", \
                      len( listLA ), len( listLB )
            else :
                for iLayer in range( len( listLA ) ) :
                    name = " L[" + str( iLayer ) + "]"
                    resultComparison = compareTwo( listLA[iLayer][0], listLA[iLayer][1],
                                                   listLB[iLayer][0], listLB[iLayer][1] )
                    if ( isLayerInformationOn ) :
                        listResults.append( ( name, listLA[iLayer][0], listLB[iLayer][0],
                                              resultComparison ) )
        elif ( i == 4 ) :
            name = " fL1"
            resultComparison = compareTwo( fL1A[0], fL1A[1], fL1B[0], fL1B[1] )
            listResults.append( ( name, fL1A[0], fL1B[0], resultComparison ) )
        elif ( i == 5 ) :
            name = " fL2"
            resultComparison = compareTwo( fL2A[0], fL2A[1], fL2B[0], fL2B[1] )
            listResults.append( ( name, fL2A[0], fL2B[0], resultComparison ) )
        elif ( i == 6 ) :
            name = " fL3"
            resultComparison = compareTwo( fL3A[0], fL3A[1], fL3B[0], fL3B[1] )
            listResults.append( ( name, fL3A[0], fL3B[0], resultComparison ) )
        elif ( i == 7 ) :
            name = " fL4"
            resultComparison = compareTwo( fL4A[0], fL4A[1], fL4B[0], fL4B[1] )
            listResults.append( ( name, fL4A[0], fL4B[0], resultComparison ) )
        elif ( i == 8 ) :
            if ( len( listRA ) != len( listRB ) ) :
                print " ***DIFFERENT NUMBER OF RINGS *** : ", \
                      len( listRA ), len( listRB )
            else :
                for iRing in range( len( listRA ) ) :
                    name = " R[" + str( iRing ) + "]"
                    resultComparison = compareTwo( listRA[iRing][0], listRA[iRing][1],
                                                   listRB[iRing][0], listRB[iRing][1] )
                    if ( isRingInformationOn ) :
                        listResults.append( ( name, listRA[iRing][0], listRB[iRing][0],
                                              resultComparison ) )
        elif ( i == 9 ) :
            name = " fR1"
            resultComparison = compareTwo( fR1A[0], fR1A[1], fR1B[0], fR1B[1] )
            listResults.append( ( name, fR1A[0], fR1B[0], resultComparison ) )
        elif ( i == 10 ) :
            name = " fR2"
            resultComparison = compareTwo( fR2A[0], fR2A[1], fR2B[0], fR2B[1] )
            listResults.append( ( name, fR2A[0], fR2B[0], resultComparison ) )
        elif ( i == 11 ) :
            name = " fR3"
            resultComparison = compareTwo( fR3A[0], fR3A[1], fR3B[0], fR3B[1] )
            listResults.append( ( name, fR3A[0], fR3B[0], resultComparison ) )
        elif ( i == 12 ) :
            name = " Evis e-/+ [MeV]"
            resultComparison = compareTwo( electronEvisA[0], electronEvisA[1],
                                           electronEvisB[0], electronEvisB[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, electronEvisA[0], electronEvisB[0],
                                      resultComparison ) )
        elif ( i == 13 ) :
            name = " Etot e-/+ [MeV]"
            resultComparison = compareTwo( electronEtotA[0], electronEtotA[1],
                                           electronEtotB[0], electronEtotB[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, electronEtotA[0], electronEtotB[0],
                                      resultComparison ) )
        elif ( i == 14 ) :
            name = " fL1 e-/e+"
            resultComparison = compareTwo( electronfL1A[0], electronfL1A[1],
                                           electronfL1B[0], electronfL1B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, electronfL1A[0], electronfL1B[0],
                                      resultComparison ) )
        elif ( i == 15 ) :
            name = " fL2 e-/e+"
            resultComparison = compareTwo( electronfL2A[0], electronfL2A[1],
                                           electronfL2B[0], electronfL2B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, electronfL2A[0], electronfL2B[0],
                                      resultComparison ) )
        elif ( i == 16 ) :
            name = " fL3 e-/e+"
            resultComparison = compareTwo( electronfL3A[0], electronfL3A[1],
                                           electronfL3B[0], electronfL3B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, electronfL3A[0], electronfL3B[0],
                                      resultComparison ) )
        elif ( i == 17 ) :
            name = " fL4 e-/e+"
            resultComparison = compareTwo( electronfL4A[0], electronfL4A[1],
                                           electronfL4B[0], electronfL4B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, electronfL4A[0], electronfL4B[0],
                                      resultComparison ) )
        elif ( i == 18 ) :
            name = " fR1 e-/e+"
            resultComparison = compareTwo( electronfR1A[0], electronfR1A[1],
                                           electronfR1B[0], electronfR1B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, electronfR1A[0], electronfR1B[0],
                                      resultComparison ) )
        elif ( i == 19 ) :
            name = " fR2 e-/e+"
            resultComparison = compareTwo( electronfR2A[0], electronfR2A[1],
                                           electronfR2B[0], electronfR2B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, electronfR2A[0], electronfR2B[0],
                                      resultComparison ) )
        elif ( i == 20 ) :
            name = " fR3 e-/e+"
            resultComparison = compareTwo( electronfR3A[0], electronfR3A[1],
                                           electronfR3B[0], electronfR3B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, electronfR3A[0], electronfR3B[0],
                                      resultComparison ) )
        elif ( i == 21 ) :
            name = " Evis protons [MeV]"
            resultComparison = compareTwo( protonEvisA[0], protonEvisA[1],
                                           protonEvisB[0], protonEvisB[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, protonEvisA[0], protonEvisB[0],
                                      resultComparison ) )
        elif ( i == 22 ) :
            name = " Etot protons [MeV]"
            resultComparison = compareTwo( protonEtotA[0], protonEtotA[1],
                                           protonEtotB[0], protonEtotB[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, protonEtotA[0], protonEtotB[0],
                                      resultComparison ) )
        elif ( i == 23 ) :
            name = " fL1 protons"
            resultComparison = compareTwo( protonfL1A[0], protonfL1A[1],
                                           protonfL1B[0], protonfL1B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, protonfL1A[0], protonfL1B[0],
                                      resultComparison ) )
        elif ( i == 24 ) :
            name = " fL2 protons"
            resultComparison = compareTwo( protonfL2A[0], protonfL2A[1],
                                           protonfL2B[0], protonfL2B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, protonfL2A[0], protonfL2B[0],
                                      resultComparison ) )
        elif ( i == 25 ) :
            name = " fL3 protons"
            resultComparison = compareTwo( protonfL3A[0], protonfL3A[1],
                                           protonfL3B[0], protonfL3B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, protonfL3A[0], protonfL3B[0],
                                      resultComparison ) )
        elif ( i == 26 ) :
            name = " fL4 protons"
            resultComparison = compareTwo( protonfL4A[0], protonfL4A[1],
                                           protonfL4B[0], protonfL4B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, protonfL4A[0], protonfL4B[0],
                                      resultComparison ) )
        elif ( i == 27 ) :
            name = " fR1 protons"
            resultComparison = compareTwo( protonfR1A[0], protonfR1A[1],
                                           protonfR1B[0], protonfR1B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, protonfR1A[0], protonfR1B[0],
                                      resultComparison ) )
        elif ( i == 28 ) :
            name = " fR2 protons"
            resultComparison = compareTwo( protonfR2A[0], protonfR2A[1],
                                           protonfR2B[0], protonfR2B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, protonfR2A[0], protonfR2B[0],
                                      resultComparison ) )
        elif ( i == 29 ) :
            name = " fR3 protons"
            resultComparison = compareTwo( protonfR3A[0], protonfR3A[1],
                                           protonfR3B[0], protonfR3B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, protonfR3A[0], protonfR3B[0],
                                      resultComparison ) )
        elif ( i == 30 ) :
            name = " Evis pi+/- [MeV]"
            resultComparison = compareTwo( pionEvisA[0], pionEvisA[1],
                                           pionEvisB[0], pionEvisB[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pionEvisA[0], pionEvisB[0],
                                      resultComparison ) )
        elif ( i == 31 ) :
            name = " Etot pi+/- [MeV]"
            resultComparison = compareTwo( pionEtotA[0], pionEtotA[1],
                                           pionEtotB[0], pionEtotB[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pionEtotA[0], pionEtotB[0],
                                      resultComparison ) )
        elif ( i == 32 ) :
            name = " fL1 pi+/-"
            resultComparison = compareTwo( pionfL1A[0], pionfL1A[1],
                                           pionfL1B[0], pionfL1B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pionfL1A[0], pionfL1B[0],
                                      resultComparison ) )
        elif ( i == 33 ) :
            name = " fL2 pi+/-"
            resultComparison = compareTwo( pionfL2A[0], pionfL2A[1],
                                           pionfL2B[0], pionfL2B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pionfL2A[0], pionfL2B[0],
                                      resultComparison ) )
        elif ( i == 34 ) :
            name = " fL3 pi+/-"
            resultComparison = compareTwo( pionfL3A[0], pionfL3A[1],
                                           pionfL3B[0], pionfL3B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pionfL3A[0], pionfL3B[0],
                                      resultComparison ) )
        elif ( i == 35 ) :
            name = " fL4 pi+/-"
            resultComparison = compareTwo( pionfL4A[0], pionfL4A[1],
                                           pionfL4B[0], pionfL4B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pionfL4A[0], pionfL4B[0],
                                      resultComparison ) )
        elif ( i == 36 ) :
            name = " fR1 pi+/-"
            resultComparison = compareTwo( pionfR1A[0], pionfR1A[1],
                                           pionfR1B[0], pionfR1B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pionfR1A[0], pionfR1B[0],
                                      resultComparison ) )
        elif ( i == 37 ) :
            name = " fR2 pi+/-"
            resultComparison = compareTwo( pionfR2A[0], pionfR2A[1],
                                           pionfR2B[0], pionfR2B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pionfR2A[0], pionfR2B[0],
                                      resultComparison ) )
        elif ( i == 38 ) :
            name = " fR3 pi+/-"
            resultComparison = compareTwo( pionfR3A[0], pionfR3A[1],
                                           pionfR3B[0], pionfR3B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pionfR3A[0], pionfR3B[0],
                                      resultComparison ) )
        elif ( i == 39 ) :
            name = " Evis nucleus [MeV]"
            resultComparison = compareTwo( nucleusEvisA[0], nucleusEvisA[1],
                                           nucleusEvisB[0], nucleusEvisB[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, nucleusEvisA[0], nucleusEvisB[0],
                                      resultComparison ) )
        elif ( i == 40 ) :
            name = " Etot nucleus [MeV]"
            resultComparison = compareTwo( nucleusEtotA[0], nucleusEtotA[1],
                                           nucleusEtotB[0], nucleusEtotB[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, nucleusEtotA[0], nucleusEtotB[0],
                                      resultComparison ) )
        elif ( i == 41 ) :
            name = " fL1 nucleus"
            resultComparison = compareTwo( nucleusfL1A[0], nucleusfL1A[1],
                                           nucleusfL1B[0], nucleusfL1B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, nucleusfL1A[0], nucleusfL1B[0],
                                      resultComparison ) )
        elif ( i == 42 ) :
            name = " fL2 nucleus"
            resultComparison = compareTwo( nucleusfL2A[0], nucleusfL2A[1],
                                           nucleusfL2B[0], nucleusfL2B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, nucleusfL2A[0], nucleusfL2B[0],
                                      resultComparison ) )
        elif ( i == 43 ) :
            name = " fL3 nucleus"
            resultComparison = compareTwo( nucleusfL3A[0], nucleusfL3A[1],
                                           nucleusfL3B[0], nucleusfL3B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, nucleusfL3A[0], nucleusfL3B[0],
                                      resultComparison ) )
        elif ( i == 44 ) :
            name = " fL4 nucleus"
            resultComparison = compareTwo( nucleusfL4A[0], nucleusfL4A[1],
                                           nucleusfL4B[0], nucleusfL4B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, nucleusfL4A[0], nucleusfL4B[0],
                                      resultComparison ) )
        elif ( i == 45 ) :
            name = " fR1 nucleus"
            resultComparison = compareTwo( nucleusfR1A[0], nucleusfR1A[1],
                                           nucleusfR1B[0], nucleusfR1B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, nucleusfR1A[0], nucleusfR1B[0],
                                      resultComparison ) )
        elif ( i == 46 ) :
            name = " fR2 nucleus"
            resultComparison = compareTwo( nucleusfR2A[0], nucleusfR2A[1],
                                           nucleusfR2B[0], nucleusfR2B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, nucleusfR2A[0], nucleusfR2B[0],
                                      resultComparison ) )
        elif ( i == 47 ) :
            name = " fR3 nucleus"
            resultComparison = compareTwo( nucleusfR3A[0], nucleusfR3A[1],
                                           nucleusfR3B[0], nucleusfR3B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, nucleusfR3A[0], nucleusfR3B[0],
                                      resultComparison ) )
        elif ( i == 48 ) :
            name = " Evis mu-/+ [MeV]"
            resultComparison = compareTwo( muonEvisA[0], muonEvisA[1],
                                           muonEvisB[0], muonEvisB[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, muonEvisA[0], muonEvisB[0],
                                      resultComparison ) )
        elif ( i == 49 ) :
            name = " Etot mu-/+ [MeV]"
            resultComparison = compareTwo( muonEtotA[0], muonEtotA[1],
                                           muonEtotB[0], muonEtotB[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, muonEtotA[0], muonEtotB[0],
                                      resultComparison ) )
        elif ( i == 50 ) :
            name = " Evis K+/- [MeV]"
            resultComparison = compareTwo( kaonEvisA[0], kaonEvisA[1],
                                           kaonEvisB[0], kaonEvisB[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, kaonEvisA[0], kaonEvisB[0],
                                      resultComparison ) )
        elif ( i == 51 ) :
            name = " Etot K+/- [MeV]"
            resultComparison = compareTwo( kaonEtotA[0], kaonEtotA[1],
                                           kaonEtotB[0], kaonEtotB[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, kaonEtotA[0], kaonEtotB[0],
                                      resultComparison ) )
        elif ( i == 52 ) :
            if ( len( listStepsA ) != len( listStepsB ) ) :
                print " ***DIFFERENT NUMBER OF STEPS CASES *** : ", \
                      len( listStepsA ), len( listStepsB )
            else :
                for iStep in range( len( listStepsA ) ) :
                    name = " #steps " + findWhichCase( iStep )    
                    resultComparison = compareTwo( listStepsA[iStep][0],
                                                   listStepsA[iStep][1],
                                                   listStepsB[iStep][0],
                                                   listStepsB[iStep][1] )
                    if ( isNumberOfStepsInformationOn ) :
                        listResults.append( ( name, listStepsA[iStep][0],
                                              listStepsB[iStep][0], resultComparison ) )
        elif ( i == 53 ) :
            if ( len( listTracksA ) != len( listTracksB ) ) :
                print " ***DIFFERENT NUMBER OF TRACKS CASES *** : ", \
                      len( listTracksA ), len( listTracksB )
            else :
                for iTrack in range( len( listTracksA ) ) :
                    name = " #tracks " + findWhichCase( iTrack )    
                    resultComparison = compareTwo( listTracksA[iTrack][0],
                                                   listTracksA[iTrack][1],
                                                   listTracksB[iTrack][0],
                                                   listTracksB[iTrack][1] )
                    if ( isNumberOfTracksInformationOn ) :
                        listResults.append( ( name, listTracksA[iTrack][0],
                                              listTracksB[iTrack][0], resultComparison ) )
        elif ( i == 54 ) :
            name = " track length e-/+ [mm]"
            resultComparison = compareTwo( electronLengthA[0], electronLengthA[1],
                                           electronLengthB[0], electronLengthB[1] )
            if ( isTrackLengthInformationOn ) :
                listResults.append( ( name, electronLengthA[0], electronLengthB[0],
                                      resultComparison ) )
        elif ( i == 55 ) :
            name = " track length pi+/- [mm]"
            resultComparison = compareTwo( pionLengthA[0], pionLengthA[1],
                                           pionLengthB[0], pionLengthB[1] )
            if ( isTrackLengthInformationOn ) :
                listResults.append( ( name, pionLengthA[0], pionLengthB[0],
                                      resultComparison ) )
        elif ( i == 56 ) :
            name = " track length protons [mm]"
            resultComparison = compareTwo( protonLengthA[0], protonLengthA[1],
                                           protonLengthB[0], protonLengthB[1] )
            if ( isTrackLengthInformationOn ) :
                listResults.append( ( name, protonLengthA[0], protonLengthB[0],
                                      resultComparison ) )
        elif ( i == 57 ) :
            name = " track length gammas [mm]"
            resultComparison = compareTwo( gammaLengthA[0], gammaLengthA[1],
                                           gammaLengthB[0], gammaLengthB[1] )
            if ( isTrackLengthInformationOn ) :
                listResults.append( ( name, gammaLengthA[0], gammaLengthB[0],
                                      resultComparison ) )
        elif ( i == 58 ) :
            name = " track length neutrons [mm]"
            resultComparison = compareTwo( neutronLengthA[0], neutronLengthA[1],
                                           neutronLengthB[0], neutronLengthB[1] )
            if ( isTrackLengthInformationOn ) :
                listResults.append( ( name, neutronLengthA[0], neutronLengthB[0],
                                      resultComparison ) )
        elif ( i == 59 ) :
            name = " exiting Ekin [MeV]"
            resultComparison = compareTwo( exitKinA[0], exitKinA[1],
                                           exitKinB[0], exitKinB[1] )
            if ( isExitingInformationOn ) :
                listResults.append( ( name, exitKinA[0], exitKinB[0], resultComparison ) )
        elif ( i == 60 ) :
            name = " exiting frac. gammas"
            resultComparison = compareTwo( exitFracGammasA[0], exitFracGammasA[1],
                                           exitFracGammasB[0], exitFracGammasB[1] )
            if ( isExitingInformationOn ) :
                listResults.append( ( name, exitFracGammasA[0], exitFracGammasB[0],
                                      resultComparison ) )
        elif ( i == 61 ) :
            name = " exiting frac. neutrons"
            resultComparison = compareTwo( exitFracNeutronsA[0], exitFracNeutronsA[1],
                                           exitFracNeutronsB[0], exitFracNeutronsB[1] )
            if ( isExitingInformationOn ) :
                listResults.append( ( name, exitFracNeutronsA[0], exitFracNeutronsB[0],
                                      resultComparison ) )
        elif ( i == 62 ) :
            name = " exiting frac. neutrinos"
            resultComparison = compareTwo( exitFracNeutrinosA[0], exitFracNeutrinosA[1],
                                           exitFracNeutrinosB[0], exitFracNeutrinosB[1] )
            if ( isExitingInformationOn ) :
                listResults.append( ( name, exitFracNeutrinosA[0], exitFracNeutrinosB[0],
                                      resultComparison ) )
        elif ( i == 63 ) :
            name = " exiting frac. mu-/+"
            resultComparison = compareTwo( exitFracMuonsA[0], exitFracMuonsA[1],
                                           exitFracMuonsB[0], exitFracMuonsB[1] )
            if ( isExitingInformationOn ) :
                listResults.append( ( name, exitFracMuonsA[0], exitFracMuonsB[0],
                                      resultComparison ) )
        elif ( i == 64 ) :
            name = " exiting frac. e-/+"
            resultComparison = compareTwo( exitFracElectronsA[0], exitFracElectronsA[1],
                                           exitFracElectronsB[0], exitFracElectronsB[1] )
            if ( isExitingInformationOn ) :
                listResults.append( ( name, exitFracElectronsA[0], exitFracElectronsB[0],
                                      resultComparison ) )
        elif ( i == 65 ) :
            name = " exitFracOthers"
            resultComparison = compareTwo( exitFracOthersA[0], exitFracOthersA[1],
                                           exitFracOthersB[0], exitFracOthersB[1] )
            if ( isExitingInformationOn ) :
                listResults.append( ( name, exitFracOthersA[0], exitFracOthersB[0],
                                      resultComparison ) )
        elif ( i == 66 ) :
            name = " #exiting"
            resultComparison = compareTwo( exitNumA[0], exitNumA[1],
                                           exitNumB[0], exitNumB[1] )
            if ( isExitingInformationOn ) :
                listResults.append( ( name, exitNumA[0], exitNumB[0], resultComparison ) )
        elif ( i == 67 ) :
            name = " #exiting gammas"
            resultComparison = compareTwo( exitNumGammasA[0], exitNumGammasA[1],
                                           exitNumGammasB[0], exitNumGammasB[1] )
            if ( isExitingInformationOn ) :
                listResults.append( ( name, exitNumGammasA[0], exitNumGammasB[0],
                                      resultComparison ) )
        elif ( i == 68 ) :
            name = " #exiting neutrons"
            resultComparison = compareTwo( exitNumNeutronsA[0], exitNumNeutronsA[1],
                                           exitNumNeutronsB[0], exitNumNeutronsB[1] )
            if ( isExitingInformationOn ) :
                listResults.append( ( name, exitNumNeutronsA[0], exitNumNeutronsB[0],
                                      resultComparison ) )
        elif ( i == 69 ) :
            name = " #exiting neutrinos"
            resultComparison = compareTwo( exitNumNeutrinosA[0], exitNumNeutrinosA[1],
                                           exitNumNeutrinosB[0], exitNumNeutrinosB[1] )
            if ( isExitingInformationOn ) :
                listResults.append( ( name, exitNumNeutrinosA[0], exitNumNeutrinosB[0],
                                      resultComparison ) )
        elif ( i == 70 ) :
            name = " #exiting mu-/+"
            resultComparison = compareTwo( exitNumMuonsA[0], exitNumMuonsA[1],
                                           exitNumMuonsB[0], exitNumMuonsB[1] )
            if ( isExitingInformationOn ) :
                listResults.append( ( name, exitNumMuonsA[0], exitNumMuonsB[0],
                                      resultComparison ) )
        elif ( i == 71 ) :
            name = " #exiting e-/+"
            resultComparison = compareTwo( exitNumElectronsA[0], exitNumElectronsA[1],
                                           exitNumElectronsB[0], exitNumElectronsB[1] )
            if ( isExitingInformationOn ) :
                listResults.append( ( name, exitNumElectronsA[0], exitNumElectronsB[0],
                                      resultComparison ) )
        elif ( i == 72 ) :
            name = " #exit others"
            resultComparison = compareTwo( exitNumOthersA[0], exitNumOthersA[1],
                                           exitNumOthersB[0], exitNumOthersB[1] )
            if ( isExitingInformationOn ) :
                listResults.append( ( name, exitNumOthersA[0], exitNumOthersB[0],
                                      resultComparison ) )
        else :
            print "CASE THAT SHOULD NOT HAPPEN : i = ", i

    count = 0
    for result in listResults :
        count += 1
        if ( ( isConditionAND  and
               ( math.fabs( result[3][0] ) > thresholdSigmaDifference  and
                 math.fabs( result[3][1] ) > thresholdRelativeDifference ) )
             or
             ( not isConditionAND  and
               ( math.fabs( result[3][0] ) > thresholdSigmaDifference  or
                 math.fabs( result[3][1] ) > thresholdRelativeDifference ) ) ) :
            print result[0],
            print '%.2f' % result[1],
            print '%.2f' % result[2],
            print " diff:", 
            print '%.2f' % result[3][0],
            print "sigma , ",
            print '%.2f' % ( result[3][1]*100.0 ),
            print "%"

    return


def findWhichCase( i ) :
    # Given the case (integer) for the average number per event of
    # either steps or tracks, it returns the string with the name
    # of the particle.
    
    caseName = " "
    if ( i == 0 ) :
        caseName = "total" 
    elif ( i == 1 ) :
        caseName = "positives" 
    elif ( i == 2 ) :
        caseName = "neutrals" 
    elif ( i == 3 ) :
        caseName = "negatives" 
    elif ( i == 4 ) :
        caseName = "particles with 0 PDG code" 
    elif ( i == 5 ) :
        caseName = "particles with Unrecognized PDG code" 
    elif ( i == 6 ) :
        caseName = "electromagnetic (e+ , e- , gammas)"
    elif ( i == 7 ) :
        caseName = "electroweak (mu+, mu-, tau+, tau-, neutrinos)"
    elif ( i == 8 ) :
        caseName = "hadrons"
    elif ( i == 9 ) :
        caseName = "mesons"
    elif ( i == 10 ) :
        caseName = "baryons"
    elif ( i == 11 ) :
        caseName = "light mesons (u/ubar/d/dbar)"
    elif ( i == 12 ) :
        caseName = "light baryons (u/ubar/d/dbar)"
    elif ( i == 13 ) :
        caseName = "strange (s/sbar) mesons"
    elif ( i == 14 ) :
        caseName = "strange (s/sbar) baryons"
    elif ( i == 15 ) :
        caseName = "heavy (c/cbar or b/bbar) mesons"
    elif ( i == 16 ) :
        caseName = "heavy (c/cbar or b/bbar) baryons"
    elif ( i == 17 ) :
        caseName = "electrons"
    elif ( i == 18 ) :
        caseName = "gammas"
    elif ( i == 19 ) :
        caseName = "positrons"
    elif ( i == 20 ) :
        caseName = "mu-"
    elif ( i == 21 ) :
        caseName = "mu+"
    elif ( i == 22 ) :
        caseName = "tau-"
    elif ( i == 23 ) :
        caseName = "tau+"
    elif ( i == 24 ) :
        caseName = "neutrinos"
    elif ( i == 25 ) :
        caseName = "pi+"
    elif ( i == 26 ) :
        caseName = "pi0"
    elif ( i == 27 ) :
        caseName = "pi-"
    elif ( i == 28 ) :
        caseName = "K+"
    elif ( i == 29 ) :
        caseName = "K-neutral (K0/K0bar or K0_S/K0_L)"
    elif ( i == 30 ) :
        caseName = "K-"
    elif ( i == 31 ) :
        caseName = "protons"
    elif ( i == 32 ) :
        caseName = "anti-protons"
    elif ( i == 33 ) :
        caseName = "neutrons"
    elif ( i == 34 ) :
        caseName = "anti-neutrons"
    else :
        print " UNKNOWN CASE ", i
        
    return caseName


def compareTwo( x, sigma_x, y, sigma_y ) :
    # Given two numbers, x and y, and their statistical errors,
    # sigma_x and sigma_y respectively, this function returns
    # the difference y - x in terms of its statistical error
    # (i.e. the number of sigma of deviation), and in terms
    # of the average of the two quantities, (y+x)/2, (i.e.
    # the relative difference).
    # Notice that these two differences can be positive or
    # negative, because they are calculated from (y - x).

    # Protect against extremely large values.
    if ( sigma_x > 1.0E+20 ) :
        print " ***WARNING*** crazy sigma_x = ", sigma_x
        sigma_x = 1.0E+20
    if ( sigma_y > 1.0E+20 ) :
        print " ***WARNING*** crazy sigma_y = ", sigma_y
        sigma_y = 1.0E+20

    sigmaDiff = 0.0
    if ( sigma_x > 0.0  or  sigma_y > 0.0 ) :
        sigma = math.sqrt( sigma_x*sigma_x + sigma_y*sigma_y )    
        sigmaDiff = ( y - x ) / sigma

    relativeDiff = 0.0
    if ( math.fabs( x ) > 0  or  math.fabs( y ) > 0 ) :
        relativeDiff = 2.0 * ( y - x ) / ( y + x ) 

    return ( sigmaDiff, relativeDiff )


#===============================================
#==================== MAIN ===================== 
#===============================================

# Check whether the number of arguments are correct or not.
if ( len( sys.argv ) != 3 ) :
    print " Usage:  compare2 fileA fileB "
else :
    fileA = open( sys.argv[1], "r" )
    fileB = open( sys.argv[2], "r" )

    valuesA = ( cpuTimeA, visEnergyA, totEnergyA,
                listLA, fL1A, fL2A, fL3A, fL4A,
                listRA, fR1A, fR2A, fR3A,
                electronEvisA, electronEtotA,
                electronfL1A, electronfL2A, electronfL3A, electronfL4A, 
                electronfR1A, electronfR2A, electronfR3A,
                protonEvisA, protonEtotA,
                protonfL1A, protonfL2A, protonfL3A, protonfL4A, 
                protonfR1A, protonfR2A, protonfR3A,
                pionEvisA, pionEtotA,
                pionfL1A, pionfL2A, pionfL3A, pionfL4A, 
                pionfR1A, pionfR2A, pionfR3A,
                nucleusEvisA, nucleusEtotA,
                nucleusfL1A, nucleusfL2A, nucleusfL3A, nucleusfL4A, 
                nucleusfR1A, nucleusfR2A, nucleusfR3A,
                muonEvisA, muonEtotA,
                kaonEvisA, kaonEtotA,
                listStepsA, listTracksA,
                electronLengthA, pionLengthA, protonLengthA,
                gammaLengthA, neutronLengthA,
                exitKinA, exitFracGammasA, exitFracNeutronsA, exitFracNeutrinosA,
                exitFracMuonsA, exitFracElectronsA, exitFracOthersA,
                exitNumA, exitNumGammasA, exitNumNeutronsA, exitNumNeutrinosA,
                exitNumMuonsA, exitNumElectronsA, exitNumOthersA
                ) = funExtract( fileA )

    valuesB = ( cpuTimeB, visEnergyB, totEnergyB,
                listLB, fL1B, fL2B, fL3B, fL4B,
                listRB, fR1B, fR2B, fR3B,
                electronEvisB, electronEtotB,
                electronfL1B, electronfL2B, electronfL3B, electronfL4B, 
                electronfR1B, electronfR2B, electronfR3B,
                protonEvisB, protonEtotB,
                protonfL1B, protonfL2B, protonfL3B, protonfL4B, 
                protonfR1B, protonfR2B, protonfR3B,
                pionEvisB, pionEtotB,
                pionfL1B, pionfL2B, pionfL3B, pionfL4B, 
                pionfR1B, pionfR2B, pionfR3B,
                nucleusEvisB, nucleusEtotB,
                nucleusfL1B, nucleusfL2B, nucleusfL3B, nucleusfL4B, 
                nucleusfR1B, nucleusfR2B, nucleusfR3B,
                muonEvisB, muonEtotB,
                kaonEvisB, kaonEtotB,
                listStepsB, listTracksB,
                electronLengthB, pionLengthB, protonLengthB,
                gammaLengthB, neutronLengthB,
                exitKinB, exitFracGammasB, exitFracNeutronsB, exitFracNeutrinosB,
                exitFracMuonsB, exitFracElectronsB, exitFracOthersB,
                exitNumB, exitNumGammasB, exitNumNeutronsB, exitNumNeutrinosB,
                exitNumMuonsB, exitNumElectronsB, exitNumOthersB
                ) = funExtract( fileB )

    skip = 0
    if ( cpuTimeA < 1.0E-9 ) :
        print " ***WARNING*** : 1st file WITHOUT RUN SUMMARY! "
        skip = 1
    if ( cpuTimeB < 1.0E-9 ) :
        print " ***WARNING*** : 2nd file WITHOUT RUN SUMMARY! "
        skip = 1
    if ( skip ) :
        print "      ---> Skipping comparison!"
    else :
        doStatisticalComparison( valuesA, valuesB )
        
    fileA.close()
    fileB.close()

#-------------------------------------------------------------------

