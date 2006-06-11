#!/usr/bin/python

#--------------------------------------------------------------------
# Last update: 11-Jun-2006
#
# This script should be run in the directory which has, as immediate
# subdirectories, the results of the Grid validation testing,
# which consist of tar-balls (one per subdirectory, each subdirectory
# being a job).
#
# The script does is following: it goes in each subdirectory,
# and look for the two output.log- files, and then compare all the
# observables that are at the end of these files:
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
# The script prints out those observables that differ by more
# than  5 sigma  and/or  by more than  5%
# ( NB) See  ***LOOKHERE***  for the choice of the thresholds and
#       and whether the condition is based on "and" or "or".)
#
# At the end, the script prints out also a summary results of the
# average discrepancies (if larger than the above thresholds),
# both as overall statistics and also by calorimeter, by
# beam particle type, by beam energy.
#
# This script has not input argument, and should be run as:
#         $  python compareAll.py
#
# NB) This script does not use the script  compare2.py , but
#     it is largerly based on it: see  ***DIFFERENCE***  to
#     look at the differences between them.
#
#--------------------------------------------------------------------

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
    pdg0Evis = (0.0, 0.0)
    pdg0Etot = (0.0, 0.0)
    pdg0fL1 = (0.0, 0.0)
    pdg0fL2 = (0.0, 0.0)
    pdg0fL3 = (0.0, 0.0)
    pdg0fL4 = (0.0, 0.0)
    pdg0fR1 = (0.0, 0.0)
    pdg0fR2 = (0.0, 0.0)
    pdg0fR3 = (0.0, 0.0)
    muonEvis = (0.0, 0.0)
    kaonEvis = (0,0, 0.0)
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
        statusShowerPdg0 = 0
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
                        elif ( line.find( "pdg0" ) > -1 ) :
                            statusShowerPdg0 = 1
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
                                
                    elif ( statusShowerPdg0 ) :
                        if ( line.find( "<E_vis> =" ) > -1 ) :
                            pdg0Evis = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                         float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                            #print "***DEBUG***  pdg0Evis = ", pdg0Evis
                        elif ( line.find( "<E_tot> =" ) > -1 ) :
                            pdg0Etot = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                         float( line.split( "=" )[1].split( "+/-" )[1].split( " " )[1] ) )
                            #print "***DEBUG***  pdg0Etot = ", pdg0Etot

                        if ( statusLfractions  and  not statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                pdg0fL1 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  pdg0fL1 = ", pdg0fL1
                            elif ( line.find( "2nd" ) > -1 ) :
                                pdg0fL2 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  pdg0fL2 = ", pdg0fL2
                            elif ( line.find( "3rd" ) > -1 ) :
                                pdg0fL3 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  pdg0fL3 = ", pdg0fL3
                            elif ( line.find( "4th" ) > -1 ) :
                                pdg0fL4 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  pdg0fL4 = ", pdg0fL4
                        elif ( not statusLfractions  and statusRfractions ) :
                            if ( line.find( "1st" ) > -1 ) :
                                pdg0fR1 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  pdg0fR1 = ", pdg0fR1
                            elif ( line.find( "2nd" ) > -1 ) :
                                pdg0fR2 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  pdg0fR2 = ", pdg0fR2
                            elif ( line.find( "3rd" ) > -1 ) :
                                pdg0fR3 = ( float( line.split( "=" )[1].split( "+/-" )[0] ), 
                                            float( line.split( "=" )[1].split( "+/-" )[1].split( "%" )[0] ) )
                                #print "***DEBUG***  pdg0fR3 = ", pdg0fR3

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
             pdg0Evis, pdg0Etot,
             pdg0fL1, pdg0fL2, pdg0fL3, pdg0fL4, 
             pdg0fR1, pdg0fR2, pdg0fR3,
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
      pdg0EvisA, pdg0EtotA,
      pdg0fL1A, pdg0fL2A, pdg0fL3A, pdg0fL4A, 
      pdg0fR1A, pdg0fR2A, pdg0fR3A,
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
      pdg0EvisB, pdg0EtotB,
      pdg0fL1B, pdg0fL2B, pdg0fL3B, pdg0fL4B, 
      pdg0fR1B, pdg0fR2B, pdg0fR3B,
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
            name = " Evis pdg0 [MeV]"
            resultComparison = compareTwo( pdg0EvisA[0], pdg0EvisA[1],
                                           pdg0EvisB[0], pdg0EvisB[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pdg0EvisA[0], pdg0EvisB[0],
                                      resultComparison ) )
        elif ( i == 40 ) :
            name = " Etot pdg0 [MeV]"
            resultComparison = compareTwo( pdg0EtotA[0], pdg0EtotA[1],
                                           pdg0EtotB[0], pdg0EtotB[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pdg0EtotA[0], pdg0EtotB[0],
                                      resultComparison ) )
        elif ( i == 41 ) :
            name = " fL1 pdg0"
            resultComparison = compareTwo( pdg0fL1A[0], pdg0fL1A[1],
                                           pdg0fL1B[0], pdg0fL1B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pdg0fL1A[0], pdg0fL1B[0],
                                      resultComparison ) )
        elif ( i == 42 ) :
            name = " fL2 pdg0"
            resultComparison = compareTwo( pdg0fL2A[0], pdg0fL2A[1],
                                           pdg0fL2B[0], pdg0fL2B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pdg0fL2A[0], pdg0fL2B[0],
                                      resultComparison ) )
        elif ( i == 43 ) :
            name = " fL3 pdg0"
            resultComparison = compareTwo( pdg0fL3A[0], pdg0fL3A[1],
                                           pdg0fL3B[0], pdg0fL3B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pdg0fL3A[0], pdg0fL3B[0],
                                      resultComparison ) )
        elif ( i == 44 ) :
            name = " fL4 pdg0"
            resultComparison = compareTwo( pdg0fL4A[0], pdg0fL4A[1],
                                           pdg0fL4B[0], pdg0fL4B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pdg0fL4A[0], pdg0fL4B[0],
                                      resultComparison ) )
        elif ( i == 45 ) :
            name = " fR1 pdg0"
            resultComparison = compareTwo( pdg0fR1A[0], pdg0fR1A[1],
                                           pdg0fR1B[0], pdg0fR1B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pdg0fR1A[0], pdg0fR1B[0],
                                      resultComparison ) )
        elif ( i == 46 ) :
            name = " fR2 pdg0"
            resultComparison = compareTwo( pdg0fR2A[0], pdg0fR2A[1],
                                           pdg0fR2B[0], pdg0fR2B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pdg0fR2A[0], pdg0fR2B[0],
                                      resultComparison ) )
        elif ( i == 47 ) :
            name = " fR3 pdg0"
            resultComparison = compareTwo( pdg0fR3A[0], pdg0fR3A[1],
                                           pdg0fR3B[0], pdg0fR3B[1] )
            if ( isShowerPerParticleInformationOn ) :
                listResults.append( ( name, pdg0fR3A[0], pdg0fR3B[0],
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

    return listResults         #***DIFFERENCE*** compare2.py returns nothing.


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

    sigmaDiff = 0.0
    if ( sigma_x > 0.0  or  sigma_y > 0.0 ) :
        sigma = math.sqrt( sigma_x*sigma_x + sigma_y*sigma_y )    
        sigmaDiff = ( y - x ) / sigma

    relativeDiff = 0.0
    if ( math.fabs( x ) > 0  or  math.fabs( y ) > 0 ) :
        relativeDiff = 2.0 * ( y - x ) / ( y + x ) 

    return ( sigmaDiff, relativeDiff )


#***DIFFERENCE*** compare2.py does not have this function.
def doUpdateResults( resultsToBeUpdated, resultsNew ) :
    # Given in input the current result, resultsToBeUpdated,
    # in the form:
    #   ( nameObservable, valueA, valueB, resultComparison )
    # and the new partial result, resultsNew, in the same form,
    # this function returns a new result, which is obtained by
    # setting 0.0 the valueA and valueB, and summing the two
    # resultComparison (because we want, at the end, to calculate
    # the average...).

    results = []
    if ( len( resultsToBeUpdated ) == 0 ) :
        results = resultsNew
    elif ( len( resultsToBeUpdated ) != len( resultsNew ) ) :
        print "***CANNOT UPDATE RESULTS*** : DIFFERENT LENGTHS:", \
              len( resultsToBeUpdated ), len( resultsNew )
    else :
        for i in range( len( resultsNew ) ) :
            if ( resultsToBeUpdated[i][0] != resultsNew[i][0] ) :
                print "***CANNOT UPDATE RESULTS*** : DIFFERENT NAMES:", \
                resultsToBeUpdated[i][0], resultsNew[i][0]       
            else :
                name = resultsToBeUpdated[i][0]
                sumSigmaDiff = resultsToBeUpdated[i][3][0] + \
                               resultsNew[i][3][0]
                sumRelativeDiff = resultsToBeUpdated[i][3][1] + \
                                  resultsNew[i][3][1]
                results.append( ( name, 0.0, 0.0, ( sumSigmaDiff, sumRelativeDiff ) ) )
            
    return results


#===============================================
#==================== MAIN =====================  #***DIFFERENCE*** compare2.py has 
#===============================================  # a quite different main.

# Besides to print out the results of the comparisons for
# each pair of files (exactly as we would do if we run
# the script compare2.py on these two files), we also
# want to print, at the end, some overall, summary
# information. To do so, we need the list of results,
# one for each observable, in the form:
#   ( nameObservable, valueA, valueB, resultComparison )
# and the counting of the number of times we have updated it,
# in such a way to compute the average of the resultComparison.
# At the end, only those observables whose average resultComparison
# is beyond the acceptance window (thresholdSigmaDifference and
# thresholdRelativeDifference) are printed out in the
# finaly grand summary.

countResults = 0
sumResults = []

dictCaloResults = { 'FeSci':[ 0, [] ] ,
                    'CuSci':[ 0, [] ],
                    'CuLAr':[ 0, [] ],
                    'WLAr':[ 0, [] ],
                    'PbSci':[ 0, [] ],
                    'PbLAr':[ 0, [] ],
                    'PbWO4':[ 0, [] ] }

dictBeamPartResults = { 'e-':[ 0, [] ] ,
                        'pi+':[ 0, [] ] ,
                        'pi-':[ 0, [] ] ,
                        'k+':[ 0, [] ] ,
                        'k-':[ 0, [] ] ,
                        'k0L':[ 0, [] ] ,
                        'p':[ 0, [] ] ,
                        'n':[ 0, [] ] }

dictBeamEResults = { '1GeV':[ 0, [] ] ,
                     '2GeV':[ 0, [] ] ,
                     '3GeV':[ 0, [] ] ,
                     '4GeV':[ 0, [] ] ,
                     '5GeV':[ 0, [] ] ,
                     '6GeV':[ 0, [] ] ,
                     '7GeV':[ 0, [] ] ,
                     '8GeV':[ 0, [] ] ,
                     '9GeV':[ 0, [] ] ,
                     '10GeV':[ 0, [] ] ,
                     '20GeV':[ 0, [] ] ,
                     '30GeV':[ 0, [] ] ,
                     '50GeV':[ 0, [] ] ,
                     '60GeV':[ 0, [] ] ,
                     '80GeV':[ 0, [] ] ,
                     '100GeV':[ 0, [] ] ,
                     '120GeV':[ 0, [] ] ,
                     '150GeV':[ 0, [] ] ,
                     '180GeV':[ 0, [] ] ,
                     '200GeV':[ 0, [] ] ,
                     '250GeV':[ 0, [] ] ,
                     '300GeV':[ 0, [] ] }

os.system( "ls -1F | grep / > listDir.txt" )
listDir = open( "listDir.txt", "r" )

for dir in listDir :
    #print dir
    saveDir = os.getcwd()
    os.chdir( dir.strip() )
    currentDir = os.getcwd()
    print " --- Directory --- ", currentDir

    os.system( "ls -1 > listFiles.txt" )
    listFiles = open( "listFiles.txt", "r" )
    fileA = 0
    nameFileA = ""
    fileB = 0
    nameFileB = ""
    for iFile in listFiles :
        #print " iFile=", iFile  ###DEBUG###
        if ( iFile.find( "output.log-" ) > -1 ) :
            #print " --- Look at the output file = ", iFile.strip(), " --- "
            if ( not nameFileA ) :
                nameFileA = iFile.strip()
                fileA = open( iFile.strip(), "r" )
            elif ( not nameFileB ) :
                nameFileB = iFile.strip()
                fileB = open( iFile.strip(), "r" )
            else :
                print "***STRANGE***: MORE THAN 2 output.log- FILES!"

    listFiles.close()
    os.system( "rm listFiles.txt" )
    os.chdir( saveDir)

    print " fileA =", nameFileA
    print " fileB =", nameFileB

    if ( fileA  and  fileB ) :

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
                    pdg0EvisA, pdg0EtotA,
                    pdg0fL1A, pdg0fL2A, pdg0fL3A, pdg0fL4A, 
                    pdg0fR1A, pdg0fR2A, pdg0fR3A,
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
                    pdg0EvisB, pdg0EtotB,
                    pdg0fL1B, pdg0fL2B, pdg0fL3B, pdg0fL4B, 
                    pdg0fR1B, pdg0fR2B, pdg0fR3B,
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

        listResultsThisCase = doStatisticalComparison( valuesA, valuesB )

        updatedResults = doUpdateResults( sumResults, listResultsThisCase )
        if ( len( updatedResults ) ) :
            countResults += 1
            sumResults = updatedResults
        
        for iDetector in dictCaloResults.keys() :
            if ( nameFileA.find( "-" + iDetector + "-" ) > -1 ) :
                updatedResults = doUpdateResults( dictCaloResults[ iDetector ][1],
                                                  listResultsThisCase )
                if ( len( updatedResults ) ) :
                    dictCaloResults[ iDetector ] = \
                                     ( dictCaloResults[ iDetector ][0] + 1,
                                       updatedResults )

        for iPart in dictBeamPartResults.keys() :
            if ( nameFileA.find( "-" + iPart + "-" ) > -1 ) :
                updatedResults = doUpdateResults( dictBeamPartResults[ iPart ][1],
                                                  listResultsThisCase )
                if ( len( updatedResults ) ) :
                    dictBeamPartResults[ iPart ] = \
                                         ( dictBeamPartResults[ iPart ][0] + 1,
                                           updatedResults )

        for iEnergy in dictBeamEResults.keys() :
            if ( nameFileA.find( "-" + iEnergy + "-" ) > -1 ) :
                updatedResults = doUpdateResults( dictBeamEResults[ iEnergy ][1],
                                                  listResultsThisCase )
                if ( len( updatedResults ) ) :
                    dictBeamEResults[ iEnergy ] = \
                                      ( dictBeamEResults[ iEnergy ][0] + 1,
                                        updatedResults )
                    
    fileA.close()
    fileB.close()

print " "
print " === SUMMARY === "
print " "
print " --- Overall --- "
print " counts =", countResults
if ( countResults ) :
    for result in sumResults :
        sigmaDifference = result[3][0] / float( countResults )
        relativeDifference = result[3][1] / float( countResults )
        if ( math.fabs( sigmaDifference ) > thresholdSigmaDifference  and
             math.fabs( relativeDifference ) > thresholdRelativeDifference ) :
            print "\t", result[0], " diff:",
            print '%.2f' % sigmaDifference,
            print "sigma , ",
            print '%.2f' % ( relativeDifference*100.0 ),
            print "%"
print " "
print " --- Calorimeter --- "
for iCalo in dictCaloResults.items() :
    print "\t --- ", iCalo[0], "---"
    counts = iCalo[1][0]
    print "\t counts =", counts
    if ( counts ) :
        for result in iCalo[1][1] :
            sigmaDifference = result[3][0] / float( counts )
            relativeDifference = result[3][1] / float( counts )    
            if ( math.fabs( sigmaDifference ) > thresholdSigmaDifference  and
                 math.fabs( relativeDifference ) > thresholdRelativeDifference ) :
                print "\t", result[0], " diff:",
                print '%.2f' % sigmaDifference,
                print "sigma , ",
                print '%.2f' % ( relativeDifference*100.0 ),
                print "%"
print " "
print " --- Beam Particles --- "
for iPart in dictBeamPartResults.items() :
    print "\t --- ", iPart[0], "---"
    counts = iPart[1][0]
    print "\t counts =", counts
    if ( counts ) :
        for result in iPart[1][1] :
            sigmaDifference = result[3][0] / float( counts )
            relativeDifference = result[3][1] / float( counts )    
            if ( math.fabs( sigmaDifference ) > thresholdSigmaDifference  and
                 math.fabs( relativeDifference ) > thresholdRelativeDifference ) :
                print "\t", result[0], " diff:",
                print '%.2f' % sigmaDifference,
                print "sigma , ",
                print '%.2f' % ( relativeDifference*100.0 ),
                print "%"
print " "
print " --- Beam Energies --- "
for iEnergy in dictBeamEResults.items() :
    print "\t --- ", iEnergy[0], "---"
    counts = iEnergy[1][0]
    print "\t counts =", counts
    if ( counts ) :
        for result in iPart[1][1] :
            sigmaDifference = result[3][0] / float( counts )
            relativeDifference = result[3][1] / float( counts )    
            if ( math.fabs( sigmaDifference ) > thresholdSigmaDifference  and
                 math.fabs( relativeDifference ) > thresholdRelativeDifference ) :
                print "\t", result[0], " diff:",
                print '%.2f' % sigmaDifference,
                print "sigma , ",
                print '%.2f' % ( relativeDifference*100.0 ),
                print "%"
print " "

listDir.close()

#-------------------------------------------------------------------

