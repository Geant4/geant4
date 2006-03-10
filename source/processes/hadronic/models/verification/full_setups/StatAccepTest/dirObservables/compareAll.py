#!/usr/bin/python

#--------------------------------------------------------------------
# Last update: 10-Mar-2006
#
# This script should be run in the directory which has, as immediate
# subdirectories, the results of the Grid validation testing,
# which consist of tar-balls (one per subdirectory, each subdirectory
# being a job).
#
# The script does is following: it goes in each subdirectory,
# and look for the two output.log- files, and then compare all the
# observables that at the end of these files:
#
#   o  CPU times (in second, without error)
#   o  total visible energy (in MeV, with error)
#   o  total energy in the whole calorimeter (in MeV, with error) 
#   o  visible energy in each layer (in MeV, with error)
#   o  visible energy in each ring (in MeV, with error)
#   o  all the information regarding the average number of steps
#      per event (adimensional number, with error)
#   o  all the information regarding the average number of tracks
#      per event (adimensional number, with error)
# 
# The script prints out those observables that diffear by more
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
#--------------------------------------------------------------------

import os
import sys
import string
import math

#***LOOKHERE***
thresholdSigmaDifference = 5.0
thresholdRelativeDifference = 0.05
isConditionAND = 0    # 1 = AND  ; 0 = OR


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

    if ( theFile ) :
        statusEndRun = 0
        statusSteps = 0
        statusTracks = 0
        for line in theFile :
            #print line
            if ( line.find( "Run Summary" ) > -1 ) :
                statusEndRun = 1

            if ( statusEndRun ) :
                if ( line.find( "User=" ) > -1 ) :
                    cpuTime = float( line.split( "s " )[0].split( "=" )[1] )
                elif ( line.find( "active layers =" ) > -1 ) :
                    visEnergy = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                  float( line.split( "=" )[1].split( "+/-" )[1] ) )
                elif ( line.find( "whole calorimeter =" ) > -1 ) :
                    totEnergy = ( float( line.split( "=" )[1].split( "+/-" )[0] ) ,
                                  float( line.split( "=" )[1].split( "+/-" )[1] ) )
                elif ( line.find( "layer =" ) > -1 ) :
                    listL.append( ( float( line.split( "=" )[2].split( "+/-" )[0] ),
                                    float( line.split( "=" )[2].split( "+/-" )[1] ) ) )
                elif ( line.find( "iBinR =" ) > -1 ) :
                    listR.append( ( float( line.split( "=" )[2].split( "+/-" )[0] ),
                                    float( line.split( "=" )[2].split( "+/-" )[1] ) ) )
                        
                if ( statusSteps ) :
                    if ( line.find( "+/-" ) > -1 ) :
                        listSteps.append( ( float( line.split( "=" )[1].split( "+/-" )[0] ) , float( line.split( "=" )[1].split( "+/-" )[1] ) ) )
                    else :
                        statusSteps = 0

                if ( statusTracks ) :
                    if ( line.find( "+/-" ) > -1 ) :
                        listTracks.append( ( float( line.split( "=" )[1].split( "+/-" )[0] ) , float( line.split( "=" )[1].split( "+/-" )[1] ) ) )
                    else :
                        statusTracks = 0

                if ( line.find( "STEPS" ) > -1 ) :
                    statusSteps = 1
                elif ( line.find( "TRACKS" ) > -1 ) :
                    statusTracks = 1

    return ( cpuTime, visEnergy, totEnergy, listL, listR, listSteps, listTracks )


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

    ( cpuTimeA, visEnergyA, totEnergyA, listLA, listRA,
      listStepsA, listTracksA ) = valuesA
    ( cpuTimeB, visEnergyB, totEnergyB, listLB, listRB,
      listStepsB, listTracksB ) = valuesB

    listResults = []  # A list of results, where the result is a tuple:
                      # (nameObservable, valueA, valueB, resultComparison)
                      # where  resultComparison  is, in turn, a tuple:
                      #   ( sigmaDiff, relativeDiff )
                      
    for i in range(7) :
        if ( i == 0 ) :
            name = " CPU times"
            resultComparison = compareTwo( cpuTimeA, 0.0, cpuTimeB, 0.0 )
            listResults.append( ( name, cpuTimeA, cpuTimeB, resultComparison ) )
        elif ( i == 1 ) :
            name = " Visible Energy"
            resultComparison = compareTwo( visEnergyA[0], visEnergyA[1],
                                           visEnergyB[0], visEnergyB[1] )
            listResults.append( ( name, visEnergyA[0], visEnergyB[0], resultComparison ) )
        elif ( i == 2 ) :
            name = " Total Energy"
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
                    listResults.append( ( name, listLA[iLayer][0], listLB[iLayer][0],
                                          resultComparison ) )
        elif ( i == 4 ) :
            if ( len( listRA ) != len( listRB ) ) :
                print " ***DIFFERENT NUMBER OF RINGS *** : ", \
                      len( listRA ), len( listRB )
            else :
                for iRing in range( len( listRA ) ) :
                    name = " R[" + str( iRing ) + "]"
                    resultComparison = compareTwo( listRA[iRing][0], listRA[iRing][1],
                                                   listRB[iRing][0], listRB[iRing][1] )
                    listResults.append( ( name, listRA[iRing][0], listRB[iRing][0],
                                          resultComparison ) )
        elif ( i == 5 ) :
            if ( len( listStepsA ) != len( listStepsB ) ) :
                print " ***DIFFERENT NUMBER OF STEPS CASES *** : ", \
                      len( listStepsA ), len( listStepsB )
            else :
                for iStep in range( len( listStepsA ) ) :
                    name = findWhichCase( iStep )    
                    resultComparison = compareTwo( listStepsA[iStep][0],
                                                   listStepsA[iStep][1],
                                                   listStepsB[iStep][0],
                                                   listStepsB[iStep][1] )
                    listResults.append( ( name, listStepsA[iStep][0],
                                          listStepsB[iStep][0], resultComparison ) )

        elif ( i == 6 ) :
            if ( len( listTracksA ) != len( listTracksB ) ) :
                print " ***DIFFERENT NUMBER OF TRACKS CASES *** : ", \
                      len( listTracksA ), len( listTracksB )
            else :
                for iTrack in range( len( listTracksA ) ) :
                    name = findWhichCase( iTrack )    
                    resultComparison = compareTwo( listTracksA[iTrack][0],
                                                   listTracksA[iTrack][1],
                                                   listTracksB[iTrack][0],
                                                   listTracksB[iTrack][1] )
                    listResults.append( ( name, listTracksA[iTrack][0],
                                          listTracksB[iTrack][0], resultComparison ) )
        else :
            print "CASE THAT SHOULD NOT HAPPEN : i = ", i

    count = 0
    isOnSteps = 0
    isOnTracks = 0
    for result in listResults :
        count += 1
        if ( count == 4 + len( listLA ) + len( listRA ) ) :
            isOnSteps = 1
        if ( count == 3 + len( listLA ) + len( listRA ) + len( listStepsA ) ) :
            isOnTracks = 1
        if ( ( isConditionAND  and
               ( math.fabs( result[3][0] ) > thresholdSigmaDifference  and
                 math.fabs( result[3][1] ) > thresholdRelativeDifference ) )
             or
             ( not isConditionAND  and
               ( math.fabs( result[3][0] ) > thresholdSigmaDifference  or
                 math.fabs( result[3][1] ) > thresholdRelativeDifference ) ) ) :
            if ( isOnSteps ) :
                print " STEPS "
                isOnSteps = 0
            elif ( isOnTracks ) :
                print " TRACKS "
                isOnTracks = 0
            print result[0], result[1], result[2], " diff:", \
                  result[3][0], "sigma , ", result[3][1]*100.0, "%"

    return listResults


def findWhichCase( i ) :
    # Given the case (integer) for the average number per event of
    # either steps or tracks, it returns the string with the name
    # of the observable.
    
    caseName = " "
    if ( i == 0 ) :
        caseName = "# total" 
    elif ( i == 1 ) :
        caseName = "# positives" 
    elif ( i == 2 ) :
        caseName = "# neutrals" 
    elif ( i == 3 ) :
        caseName = "# negatives" 
    elif ( i == 4 ) :
        caseName = "# particles with 0 PDG code" 
    elif ( i == 5 ) :
        caseName = "# particles with Unrecognized PDG code" 
    elif ( i == 6 ) :
        caseName = "# electromagnetic (e+ , e- , gammas)"
    elif ( i == 7 ) :
        caseName = "# electroweak (mu+, mu-, tau+, tau-, neutrinos)"
    elif ( i == 8 ) :
        caseName = "# hadrons"
    elif ( i == 9 ) :
        caseName = "# mesons"
    elif ( i == 10 ) :
        caseName = "# baryons"
    elif ( i == 11 ) :
        caseName = "# light mesons (u/ubar/d/dbar)"
    elif ( i == 12 ) :
        caseName = "# light baryons (u/ubar/d/dbar)"
    elif ( i == 13 ) :
        caseName = "# strange (s/sbar) mesons"
    elif ( i == 14 ) :
        caseName = "# strange (s/sbar) baryons"
    elif ( i == 15 ) :
        caseName = "# heavy (c/cbar or b/bbar) mesons"
    elif ( i == 16 ) :
        caseName = "# heavy (c/cbar or b/bbar) baryons"
    elif ( i == 17 ) :
        caseName = "# electrons"
    elif ( i == 18 ) :
        caseName = "# gammas"
    elif ( i == 19 ) :
        caseName = "# positrons"
    elif ( i == 20 ) :
        caseName = "# mu-"
    elif ( i == 21 ) :
        caseName = "# mu+"
    elif ( i == 22 ) :
        caseName = "# tau-"
    elif ( i == 23 ) :
        caseName = "# tau+"
    elif ( i == 24 ) :
        caseName = "# neutrinos"
    elif ( i == 25 ) :
        caseName = "# pi+"
    elif ( i == 26 ) :
        caseName = "# pi0"
    elif ( i == 27 ) :
        caseName = "# pi-"
    elif ( i == 28 ) :
        caseName = "# K+"
    elif ( i == 29 ) :
        caseName = "# K-neutral (K0/K0bar or K0_S/K0_L)"
    elif ( i == 30 ) :
        caseName = "# K-"
    elif ( i == 31 ) :
        caseName = "# protons"
    elif ( i == 32 ) :
        caseName = "# anti-protons"
    elif ( i == 33 ) :
        caseName = "# neutrons"
    elif ( i == 34 ) :
        caseName = "# anti-neutrons"
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


def doUpdateResults( resultsToBeUpdated, resultsNew ) :
    # Given in input the current result, resultsToBeUpdated
    # in the form:
    #   ( nameObservable, valueA, valueB, resultComparison )
    # and the new partial result, resultsNew, in the same form,
    # this function returns a new result, which is obtained by
    # summing the two resultComparison.

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
#==================== MAIN ===================== 
#===============================================

# We need the list of results, one for each observable, in the form:
#   ( nameObservable, valueA, valueB, resultComparison )
# and the counting of the number of times we have updated it,
# in such a way to compute the average.

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

        valuesA = ( cpuTimeA, visEnergyA, totEnergyA, listLA, listRA,
                    listStepsA, listTracksA ) = funExtract( fileA )
        valuesB = ( cpuTimeB, visEnergyB, totEnergyB, listLB, listRB,
                    listStepsB, listTracksB ) = funExtract( fileB )

        listResultsThisCase = doStatisticalComparison( valuesA, valuesB )

        #print " cpuTimeA = ", cpuTimeA
        #print " cpuTimeB = ", cpuTimeB
        #print " visEnergyA = ", visEnergyA 
        #print " visEnergyB = ", visEnergyB
        #print " totEnergyA = ", totEnergyA 
        #print " totEnergyB = ", totEnergyB
        #print " listLA = ", listLA, " ; length = ", len( listLA )
        #print " listLB = ", listLB, " ; length = ", len( listLB )
        #print " listRA = ", listRA, " ; length = ", len( listRA )
        #print " listRB = ", listRB, " ; length = ", len( listRB )
        #print " listStepsA = ", listStepsA, " ; length = ", len( listStepsA )
        #print " listStepsB = ", listStepsB, " ; length = ", len( listStepsB )
        #print " listTracksA = ", listTracksA, " ; length = ", len( listTracksA )
        #print " listTracksB = ", listTracksB, " ; length = ", len( listTracksB )

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
            print result[0], " diff:", sigmaDifference, "sigma , ", \
                  relativeDifference*100.0, "%"
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
                print "\t", result[0], " diff:", sigmaDifference, "sigma , ", \
                      relativeDifference*100.0, "%"
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
                print "\t", result[0], " diff:", sigmaDifference, "sigma , ", \
                      relativeDifference*100.0, "%"
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
                print "\t", result[0], " diff:", sigmaDifference, "sigma , ", \
                      relativeDifference*100.0, "%"
print " "

listDir.close()

#-------------------------------------------------------------------

