#!/usr/bin/python

#----------------------------------------------------------------
# 29-Apr-2005 A.R. 
#
# This Python script is used for post-processing analysis, i.e.
# to produce plots (in PostScript format) of calorimeter
# observables, like energy resolution, sampling fraction, and
# e/pi as a function of the beam energy, and shower profiles
# (longitudinal and transverse). These plots are obtained from
# the information printed out at the end of the run of a
# simulation job like:   mainStatAccepTest  run.g4
# We are assuming here that the log files of these simulation
# runs are collected in a directory, eventually with a
# subdirectory structure.
#
# This script does not have input arguments, but you have to set
# appropriately some parameters, that you can find below (search
# for the keyword  ***LOOKHERE*** ). Shortly, these are the
# parameters you need to set:
#    1)  directory : where to look for the log files from which
#                    the information are extracted. It looks in
#                    the specified directory and, recursively,
#                    to all its subdirectories.
#    2)  tupleG4Versions : the two Geant4 versions that you
#                          want to compare. A single version
#                          is also accepted (in this case, the
#                          plots are produced for this G4 version,
#                          of course without any comparison).
#    3)  tuplePhysicsLists : Geant4 Physics Lists to be considered.
#    4)  tupleCaloTypes : calorimeter types to be considered.
#    5)  tupleParticles : beam particle types to be considered.
#    6)  tupleEnergies  : beam energies to be considered.
#    7)  tupleEvents    : number of simulated events for each Run.
#
# This script can be placed in any directory, and it needs to have
# in the same directory the following kumacs , that are invoked by
# the script:
#    -  resolution.kumac
#    -  sampling_fraction.kumac
#    -  ratio_e_pi.kumac
#    -  longitudinal_profile.kumac
#    -  transverse_profile.kumac
# Notice that these kumacs are executed in batch mode
# (paw -w 0 -b name.kumac), i.e. non-interactively, without
# opening the PAW window. For this reason you would get the
# following harmless error messages:
#    ***** ERROR in IACWK : Workstation is not open 
#    ***** ERROR in ISWKWN : Invalid workstation window parameters 
#    ***** ERROR in ISWKVP : Invalid workstation window parameters 
# Please ignore these messages!
#
# This script produces in output, in the same directory where this
# script is located, the following files:
#    o  listFoundFiles.txt : list of expected simulation log files
#                            that have been found.
#    o  listMissingFiles.txt : list of expected simulation log files
#                              that have not been found.
#    o  listCreatedFiles.txt : the complete list of all the files
#                              created by this script;
#    o  listNonRetrievedFiles.txt : list of files that have been created
#                                   by this script, then closed, and then
#                                   they should have been opened again
#                                   (in reading mode), but they have not
#                                   been retrieved: of course, this should
#                                   never happens, and so this file should
#                                   always be empty if everything is ok.
#    o  The beam energy values, used by the kumacs to plot some quantities
#       as a function of the beam energy.
#    o  A set of ascii files, extracted from the list of found
#       simulation log files, whose names are:
#          -  energy_resolutions.txt-LABEL
#          -  sampling_fractions.txt-LABEL
#          -  ratios_e_pi.txt-LABEL
#          -  longitudinal_profile.txt-LABEL
#          -  transverse_profile.txt-LABEL
#       where LABEL is a string that summarizes which case
#       it corresponds (e.g.: 6.2.p02-QGSP-CuLAr-pi+-20GeV-5000 ).
#    o  A set of PostScript (.ps) files, produced by the kumacs
#       above, which plots the following observables:
#          -  resolution.ps-LABEL : energy resolution versus
#                                   beam energy;
#          -  sampling_fraction.ps-LABEL : sampling fractions
#                                   versus beam energy;
#          -  ratio_e_pi.ps-LABEL : e/pi versus beam energy;
#          -  longitudinal_profile.ps-LABEL : longitudinal
#                                   shower profile, i.e. energy
#                                   deposition versus layer number;
#          -  transverse_profile.ps-LABEL : transverse shower
#                                   profile, i.e. energy deposition
#                                   versus radius beam.
#       where LABEL is a string that summarizes which case it
#       corresponds (e.g.: 6.2.p02-7.0-QGSP-CuLAr-pi+-20GeV-5000 ).
#
#----------------------------------------------------------------

import os
import sys
import string

#***LOOKHERE***
# Look for the files in this directory and recursively in all
# its subdirectories.
#directory = "/afs/cern.ch/sw/geant4/stat_testing/g70_slc3"
directory = "."                 

# Prepare all the cases.
tupleG4Versions   = ("6.2.p02", "7.0.cand04")
#tupleG4Versions   = ("6.2.p02",)

#tuplePhysicsLists = ("LHEP", "QGSP", "QGSC", "QGSP_BIC", "QGSP_BERT")
tuplePhysicsLists = ("QGSP",)

#tupleCaloTypes    = ("FeSci", "CuSci", "PbSci", "CuLAr", "PbLAr", "WLAr", "PbWO4")
tupleCaloTypes    = ("CuLAr",)

#tupleParticles    = ("e-", "pi+", "pi-", "k+", "k-", "k0L", "p", "n")
tupleParticles    = ("e-", "pi+")

tupleEnergies     = ("1GeV", "2GeV", "3GeV", "4GeV", "5GeV", "6GeV", "7GeV",
                     "8GeV", "9GeV", "10GeV", "20GeV", "30GeV", "40GeV",
                     "50GeV", "60GeV", "80GeV", "100GeV", "120GeV", "150GeV",
#                     "180GeV", "200GeV", "250GeV", "300GeV", "1000GeV")
                     "180GeV", "200GeV", "250GeV", "300GeV")
#tupleEnergies     = ("20GeV",)

tupleEvents       = ("5000",)

#***endLOOKHERE***

# Collect in files the following lists:
#  -  the list of files that have been found;
#  -  the list of files that have not been found;
#  -  the list of files that have not been retrieved;
#  -  the list of the other files created by this script;
foundFiles = open( "listFoundFiles.txt", "w" )
missingFiles = open( "listMissingFiles.txt", "w" )
createdFiles = open( "listCreatedFiles.txt", "w" )
nonRetrievedFiles = open( "listNonRetrievedFiles.txt", "w" )

#-------------------------------------------------------------
# ---------------------   FUNCTIONS   ------------------------
#-------------------------------------------------------------

def printParameters() :
    # 
    # This function prints out the values of the various parameters,
    # and also produces the file  beam_energies.txt  which is then
    # used by the kumacs for making the plots.

    print '  --- Start function  printParameters  --- '

    print '  directory = ', directory

    print '  Geant4 versions : '
    for iG4 in tupleG4Versions :
        print '                    ', iG4

    print '  Physics Lists : '
    for iPL in tuplePhysicsLists :
        print '                  ' , iPL

    print '  Calorimeter types : '
    for iCalo in tupleCaloTypes : 
        print '                      ', iCalo

    print '  Particle types : '
    for iParticle in tupleParticles :
        print '                   ', iParticle

    fileEnergyValues = open( "beam_energies.txt", "w" )
    createdFiles.write( "beam_energies.txt" + "\n" )
    print '  Beam Energy values : '
    for iE in tupleEnergies :
        print '                       ', iE
        E = float( iE.replace("GeV","") )        
        fileEnergyValues.write( str( E ) + "\n" )
    fileEnergyValues.close()

    print '  Number of Events : '
    for iN in tupleEvents :
        print '                     ', iN

    print '  --- End   function  printParameters  --- '
    return


def extractInfo( theFile , label , fileResolution , fileSampling , fileRatio ) :
    # 
    # This function reads the file  theFile  and then extracts the
    # information useful for the various observables:
    #   -  energy resolution : the information is written in the file
    #                          fileResolution ;
    #   -  sampling fraction : the information is written in the file
    #                          fileSampling ;
    #   -  e/pi ratio : the information is written in the file  fileRatio ;
    #   -  longitudinal shower shape : a file is created and filled to store
    #                                  the visible energy for each layer;
    #   -  transverse shower shape : a file is created and filled to store
    #                                the visible energy for each radial bin.
       
    print '  --- Start function  extractInfo  --- '

    fileLongitudinal = open( "longitudinal_profile.txt-" + label , "w" )
    createdFiles.write( "longitudinal_profile.txt-" + label + "\n" )
    fileTransverse = open( "transverse_profile.txt-" + label , "w" )
    createdFiles.write( "transverse_profile.txt-" + label + "\n" )

    statusRunSummary = 0
    statusLongitudinalProfile = 0
    statusTransverseProfile = 0
    statusVisibleEnergy = 0
    for line in theFile :
        #print " line : ", line
        if ( ( not statusRunSummary )  and  line.find("Run Summary") > -1 ) :
            statusRunSummary = 1
            #print " === FOUND RUN SUMMARY === ", line #***DEBUG***
        if ( statusRunSummary  and  line.find("in each Layer") > -1 ) :
            statusLongitudinalProfile = 1
        if ( statusRunSummary  and  line.find("in each Radius bin") > -1 ) :
            statusLongitudinalProfile = 0
            statusTransverseProfile   = 1
        if ( statusRunSummary  and  line.find("Visible energy information") > -1 ) :
            statusTransverseProfile = 0
            statusVisibleEnergy     = 1
        if ( statusLongitudinalProfile  and  line.find("layer = ") > -1 ) :
            fileLongitudinal.write( line.split()[5] + "\t" + line.split()[7] + "\n" )  
            #print ' longitudinal line: ' , line.split()[5] , ' +/- ', line.split()[7]
        if ( statusTransverseProfile  and  line.find("iBinR = ") > -1 ) :
            fileTransverse.write( line.split()[5] + "\t" + line.split()[7] + "\n" )
            #print ' transverse line: ' , line.split()[5] , ' +/- ', line.split()[7]
        if ( statusVisibleEnergy  and  line.find("mu_Evis") > -1 ) :
            fileRatio.write( line.split()[2] + "\n" )
            #print ' mu_Evis line: ', line.split()[2]
        if ( statusVisibleEnergy  and  line.find("sigma_Evis") > -1 ) :
            pass;
            #print ' sigma_Evis line: ', line.split()[2]
        if ( statusVisibleEnergy  and  line.find("energy resolution") > -1 ) :
            fileResolution.write( line.split()[3] + "\n" )
            #print ' energy resolution line: ', line.split()[3]
        if ( statusVisibleEnergy  and  line.find("sampling fraction") > - 1) :
            fileSampling.write( line.split()[3] + "\n" )
            #print ' sampling fraction line: ', line.split()[3]
            statusVisibleEnergy = 0
            
    fileLongitudinal.close()
    fileTransverse.close()
    
    print '  --- End   function  extractInfo  --- '
    return


def makeShowerPlots( labelShort , g4version_a , g4version_b ) :
    # 
    # This function makes the  .ps  plots of the shower shape
    # by calling the following  kumacs :
    #   -  longitudinal_profile.kumac : energy deposit vs layer number ; 
    #   -  transverse_profile.kumac : energy deposit vs radius bin .
    # These kumacs requires in input one or two ascii files, with
    # suffix  .txt-a  and  .txt-b  . These files are symbolic links to
    # files created by the previous function  extractInfo .
    
    print '  --- Start function  makeShowerPlots  --- '
    
    #print '  labelShort  = ', labelShort
    #print '  g4version_a = ', g4version_a
    #print '  g4version_b = ', g4version_b

    # Longitudinal shower shape: energy deposit vs. layer number.
    isThereFileLongitudinal_a = 0
    isThereFileLongitudinal_b = 0
    fileName_longitudinal_a = "longitudinal_profile.txt-" + g4version_a + "-" + labelShort
    if ( os.path.isfile( fileName_longitudinal_a ) ) :
        isThereFileLongitudinal_a = 1
        commandString = "ln -sf " + fileName_longitudinal_a + " longitudinal_profile.txt-a"
        #print '  commandString-a = ', commandString
        os.system( commandString )
    else :
        nonRetrievedFiles.write( fileName_longitudinal_a + "\n" )
        print '  fileName_longitudinal_a  NOT  found! '

    fileName_longitudinal_b = "longitudinal_profile.txt-" + g4version_b + "-" + labelShort
    if ( g4version_b ) :
        commandString = "ln -sf " + fileName_longitudinal_b + " longitudinal_profile.txt-b"
        if ( os.path.isfile( fileName_longitudinal_b ) ) :
            isThereFileLongitudinal_b = 1
            #print '  commandString-b = ', commandString
            os.system( commandString )
        else :
            nonRetrievedFiles.write( fileName_longitudinal_b + "\n" )
            print '  fileName_longitudinal_b  NOT  found! '
            
    if ( isThereFileLongitudinal_a   or   isThereFileLongitudinal_b ) :
        sys.stdout.flush()   # Flush the buffer to keep the output in synch with PAW
        os.system( "paw -w 0 -b longitudinal_profile.kumac" )
        label = ""
        if ( isThereFileLongitudinal_a ) :
            label = label + g4version_a + "-"
        if ( isThereFileLongitudinal_b ) :
            label = label + g4version_b + "-"
        label = label + labelShort    
        os.system( "mv longitudinal_profile.ps longitudinal_profile.ps-" + label ) 
        createdFiles.write( "longitudinal_profile.ps-" + label + "\n" )        

    # Transverse shower shape: energy deposit vs. radius bin.
    isThereFileTransverse_a = 0
    isThereFileTransverse_b = 0
    fileName_transverse_a = "transverse_profile.txt-" + g4version_a + "-" + labelShort
    if ( os.path.isfile( fileName_transverse_a ) ) :
        isThereFileTransverse_a = 1
        commandString = "ln -sf " + fileName_transverse_a + " transverse_profile.txt-a"
        #print '  commandString-a = ', commandString
        os.system( commandString )
    else :
        nonRetrievedFiles.write( fileName_transverse_a + "\n" )
        print '  fileName_transverse_a  NOT  found! '
        
    fileName_transverse_b = "transverse_profile.txt-" + g4version_b + "-" + labelShort
    if ( g4version_b ) :
        commandString = "ln -sf " + fileName_transverse_b + " transverse_profile.txt-b"
        if ( os.path.isfile( fileName_longitudinal_b ) ) :
            isThereFileTransverse_b = 1
            #print '  commandString-b = ', commandString
            os.system( commandString )
        else :
            nonRetrievedFiles.write( fileName_transverse_b + "\n" )
            print '  fileName_transverse_b  NOT  found! '

    if ( isThereFileTransverse_a   or   isThereFileTransverse_b ) :
        sys.stdout.flush()   # Flush the buffer to keep the output in synch with PAW
        os.system( "paw -w 0 -b transverse_profile.kumac" )
        label = ""
        if ( isThereFileTransverse_a ) :
            label = label + g4version_a + "-"
        if ( isThereFileTransverse_b ) :
            label = label + g4version_b + "-"
        label = label + labelShort    
        os.system( "mv transverse_profile.ps transverse_profile.ps-" + label ) 
        createdFiles.write( "transverse_profile.ps-" + label + "\n" )        

    print '  --- End   function  makeShowerPlots  --- '
    return


def makeResolutionPlots( label , g4version_a , g4version_b ) :
    # 
    # This function makes the  .ps  plots of the energy resolution
    # and sampling fraction, by calling the following  kumacs :
    #   -  resolution.kumac : energy resolution vs. beam energy ;
    #   -  sampling_fraction.kumac : sampling fraction vs. beam energy .
    # These kumacs require in input one or two ascii files, with
    # suffix  .txt-a  and  .txt-b  . These files are symbolic links to
    # files created before by the same python script.
    
    print '  --- Start function  makeResolutionPlots  --- '
    
    #print '  label       = ', label
    #print '  g4version_a = ', g4version_a
    #print '  g4version_b = ', g4version_b

    isThereFileResolution_a = 0
    isThereFileSampling_a   = 0
    isThereFileResolution_b = 0
    isThereFileSampling_b   = 0

    theLabel = g4version_a + "-" + label
    stringFileResolution_a = "energy_resolutions.txt-" + theLabel 
    stringFileSampling_a = "sampling_fractions.txt-" + theLabel
    if ( os.path.isfile( stringFileResolution_a ) ) :
        isThereFileResolution_a = 1
        #print '  fileResolution_a  FOUND '
        commandString = "ln -sf " + stringFileResolution_a + " energy_resolutions.txt-a"
        #print '  commandString-a = ', commandString
        os.system( commandString )
    else :
        nonRetrievedFiles.write( stringFileResolution_a + "\n" )
        print '  fileResolution_a  NOT  found! '
    if ( os.path.isfile( stringFileSampling_a ) ) :        
        isThereFileSampling_a = 1
        #print '  fileSampling_a  FOUND '
        commandString = "ln -sf " + stringFileSampling_a + " sampling_fractions.txt-a"
        #print '  commandString-a = ', commandString
        os.system( commandString )
    else :
        nonRetrievedFiles.write( stringFileSampling_a + "\n" )
        print '  fileSampling_a  NOT  found!'

    stringFileResolution_b = "" 
    stringFileSampling_b = ""
    if ( g4version_b ) :
        theLabel = g4version_b + "-" + label
        stringFileResolution_b = "energy_resolutions.txt-" + theLabel 
        stringFileSampling_b = "sampling_fractions.txt-" + theLabel
        if ( os.path.isfile( stringFileResolution_b ) ) :
            isThereFileResolution_b = 1
            #print '  fileResolution_b  FOUND '
            commandString = "ln -sf " + stringFileResolution_b + " energy_resolutions.txt-b"
            #print '  commandString-b = ', commandString
            os.system( commandString )
        else :
            nonRetrievedFiles.write( stringFileResolution_b + "\n" )
            print '  fileResolution_b  NOT  found! '
        if ( os.path.isfile( stringFileSampling_b ) ) :
            isThereFileSampling_b = 1
            #print '  fileSampling_b  FOUND '
            commandString = "ln -sf " + stringFileSampling_b + " sampling_fractions.txt-b"
            #print '  commandString-b = ', commandString
            os.system( commandString )
        else :
            nonRetrievedFiles.write( stringFileSampling_b + "\n" )
            print '  fileSampling_b  NOT  found!'

    # Call the kumac to make the plot of the energy resolution vs. beam energy.
    if ( isThereFileResolution_a ) :
        sys.stdout.flush()   # Flush the buffer to keep the output in synch with PAW
        os.system( "paw -w 0 -b resolution.kumac" )
        theLabel = g4version_a + "-"
        if ( isThereFileResolution_b ) :
            theLabel = theLabel + g4version_b + "-"
        theLabel = theLabel + label
        #print '   theLabel for resolution.ps- : ', theLabel
        os.system( "mv resolution.ps resolution.ps-" + theLabel )
        createdFiles.write( "resolution.ps-" + theLabel + "\n" )        
     
    # Call the kumac to make the plot of the sampling fraction vs. beam energy.
    if ( isThereFileSampling_a ) :
        sys.stdout.flush()   # Flush the buffer to keep the output in synch with PAW
        os.system( "paw -w 0 -b sampling_fraction.kumac" )
        theLabel = g4version_a + "-"
        if ( isThereFileSampling_b ) :
            theLabel = theLabel + g4version_b + "-"
        theLabel = theLabel + label
        #print '   theLabel for sampling_fraction.ps- : ', theLabel
        os.system( "mv sampling_fraction.ps sampling_fraction.ps-" + theLabel )
        createdFiles.write( "sampling_fraction.ps-" + theLabel + "\n" )        

    print '  --- End   function  makeResolutionPlots  --- '
    return


def makeRatioPlots( label , strNumEvents , g4version_a , g4version_b ) :
    # 
    # This function makes the  .ps  plot of the e/pi ratio, by
    # calling the following  kumac :
    #   -  ratio_e_pi.kumac : e/pi vs. beam energy ;
    # This kumac requires in input one or two ascii files, with
    # suffix  .txt-a  and  .txt-b  . These files are symbolic links to
    # files created by this function, using as input several ascii files
    # created by the same python script, corresponding to the visible
    # energy for e- and pi+, for different beam energies.
    #

    print '  --- Start function  makeRatioPlots  --- '
    
    #print '  label        = ', label
    #print '  strNumEvents = ', strNumEvents
    #print '  g4version_a  = ', g4version_a
    #print '  g4version_b  = ', g4version_b

    # Open the files for e- and pi+ .
    fileRatio_e_a = 0
    fileRatio_pi_a = 0
    theLabel = g4version_a + "-" + label + "-" + "e-" + "-" + strNumEvents
    #print '  fileRatio_e_a : theLabel = ', theLabel
    if ( os.path.isfile( "ratio_e_pi.txt-" + theLabel ) ) :
        fileRatio_e_a = open( "ratio_e_pi.txt-" + theLabel , "r" )
        #print '  fileRatio_e_a  FOUND '
    else :
        nonRetrievedFiles.write( "ratio_e_pi.txt-" + theLabel + "\n" )
        print '  fileRatio_e_a  NOT  found! '
    theLabel = g4version_a + "-" + label + "-" + "pi+" + "-" + strNumEvents
    #print '  fileRatio_pi_a : theLabel = ', theLabel
    if ( os.path.isfile( "ratio_e_pi.txt-" + theLabel ) ) :
        fileRatio_pi_a = open( "ratio_e_pi.txt-" + theLabel , "r" )
        #print '  fileRatio_pi_a  FOUND '
    else :
        nonRetrievedFiles.write( "ratio_e_pi.txt-" + theLabel + "\n" )
        print '  fileRatio_pi_a  NOT  found! '

    fileRatio_e_b = 0
    fileRatio_pi_b = 0
    if ( g4version_b ) :
        theLabel = g4version_b + "-" + label + "-" + "e-" + "-" + strNumEvents
        #print '  fileRatio_e_b : theLabel = ', theLabel
        if ( os.path.isfile( "ratio_e_pi.txt-" + theLabel ) ) :
            fileRatio_e_b = open( "ratio_e_pi.txt-" + theLabel , "r" )
            #print '  fileRatio_e_b  FOUND '
        else :
            nonRetrievedFiles.write( "ratio_e_pi.txt-" + theLabel + "\n" )
            print '  fileRatio_e_b  NOT  found! '
        theLabel = g4version_b + "-" + label + "-" + "pi+" + "-" + strNumEvents
        #print '  fileRatio_pi_b : theLabel = ', theLabel
        if ( os.path.isfile( "ratio_e_pi.txt-" + theLabel ) ) :
            fileRatio_pi_b = open( "ratio_e_pi.txt-" + theLabel , "r" )
            #print '  fileRatio_pi_b  FOUND '
        else :
            nonRetrievedFiles.write( "ratio_e_pi.txt-" + theLabel + "\n" )            
            print '  fileRatio_pi_b  NOT  found! '

    # Create the file of the ratio e/pi , from the file for e- (numerator)
    # and the file for pi+ (denominator).
    isThereFileRatio_e_pi_a = 0
    stringFileRatio_e_pi_a = ""
    if ( fileRatio_e_a  and  fileRatio_pi_a ) :
        theLabel = g4version_a + "-" + label + "-" + strNumEvents
        stringFileRatio_e_pi_a = "ratio_e_pi.txt-" + theLabel 
        fileRatio_e_pi_a = open( stringFileRatio_e_pi_a , "w" )
        createdFiles.write( stringFileRatio_e_pi_a + "\n" )
        # Store the lines of the two input files, corresponding to the
        # visible energy for e- and pi+ for a given beam energy,
        # in tuples; then check that they have the same size, and then
        # fill the new file with the their ratio.
        tupleE_a = fileRatio_e_a.readlines()
        tuplePi_a = fileRatio_pi_a.readlines()
        if ( len( tupleE_a ) == len( tuplePi_a ) ) :
            isThereFileRatio_e_pi_a = 1
            for i in xrange( len( tupleE_a ) ) :
                numerator = float( tupleE_a[i].split()[0] )
                denominator = float( tuplePi_a[i].split()[0] )
                #print '  numerator=', numerator, ' denominator=', denominator
                fileRatio_e_pi_a.write( str( numerator/denominator ) + "\n" )
        else :
            print ' ERROR: E and PI files of DIFFERENT SIZE! ', \
                  len( tupleE_a ), '  ', len( tuplePi_a )
        fileRatio_e_pi_a.close()

    isThereFileRatio_e_pi_b = 0
    stringFileRatio_e_pi_b = ""
    if ( fileRatio_e_b  and  fileRatio_pi_b ) :
        theLabel = g4version_b + "-" + label + "-" + strNumEvents
        stringFileRatio_e_pi_b = "ratio_e_pi.txt-" + theLabel 
        fileRatio_e_pi_b = open( stringFileRatio_e_pi_b , "w" )
        createdFiles.write( stringFileRatio_e_pi_b + "\n" )
        tupleE_b = fileRatio_e_b.readlines()
        tuplePi_b = fileRatio_pi_b.readlines()
        if ( len( tupleE_b ) == len( tuplePi_b ) ) :
            isThereFileRatio_e_pi_b = 1
            for i in xrange( len( tupleE_b ) ) :
                numerator = float( tupleE_b[i].split()[0] )
                denominator = float( tuplePi_b[i].split()[0] )
                #print '  numerator=', numerator, ' denominator=', denominator
                fileRatio_e_pi_b.write( str( numerator/denominator ) + "\n" )
        else :
            print ' ERROR: E and PI files of DIFFERENT SIZE! ', \
                  len( tupleE_b ), '  ', len( tuplePi_b )
        fileRatio_e_pi_b.close()
        
    # Call the kumac to make the plot of the Ratio e/pi vs. beam energy.
    if ( isThereFileRatio_e_pi_a ) :
        commandString = "ln -sf " + stringFileRatio_e_pi_a + " ratios_e_pi.txt-a"
        #print '  commandString-a = ', commandString
        os.system( commandString )
        theLabel = g4version_a + "-"
        if ( isThereFileRatio_e_pi_b ) :
            commandString = "ln -sf " + stringFileRatio_e_pi_b + " ratios_e_pi.txt-b"
            #print '  commandString-b = ', commandString
            os.system( commandString )
            theLabel = theLabel + g4version_b + "-"
        theLabel = theLabel + label + "-" + strNumEvents
        #print '   theLabel for ratio_e_pi.ps- : ', theLabel
        sys.stdout.flush()   # Flush the buffer to keep the output in synch with PAW
        os.system( "paw -w 0 -b ratio_e_pi.kumac" )            
        os.system( "mv ratio_e_pi.ps ratio_e_pi.ps-" + theLabel )
        createdFiles.write( "ratio_e_pi.ps-" + theLabel + "\n" )        

    # Close the files.
    if ( fileRatio_e_a ) :
        fileRatio_e_a.close()
    if ( fileRatio_pi_a ) :
        fileRatio_pi_a.close()
    if ( fileRatio_e_b ) :
        fileRatio_e_b.close()
    if ( fileRatio_pi_b ) :
        fileRatio_pi_b.close()

    print '  --- End   function  makeRatioPlots  --- '
    return


#-------------------------------------------------------------
# ---------------------   MAIN   -----------------------------
#-------------------------------------------------------------

print '  ========== START observables.py ========== '

# Print the parameters and calculate the bins for the beam energy.
printParameters();

# Keep track of all the files created by this script.
createdFiles.write( "listFoundFiles.txt" + " \n" )
createdFiles.write( "listMissingFiles.txt" + " \n" )
createdFiles.write( "listCreatedFiles.txt" + " \n" )
createdFiles.write( "listNonRetrievedFiles.txt" + " \n" )
createdFiles.write( "thePath.log" + " \n" )

# Prepare two strings which contains the (max) two Geant4 versions.
g4version1 = tupleG4Versions[0]
g4version2 = ""
if ( len( tupleG4Versions ) > 1 ) :
    g4version2 = tupleG4Versions[1]

# Loop over all cases:
for iPL in tuplePhysicsLists :        # Loop over Physics Lists
    #print " Physics List = ", iPL

    for iCalo in tupleCaloTypes :       # Loop over Calorimeter types
        #print " Calo Type = " , iCalo

        for iN in tupleEvents :           # Loop over Event numbers
            #print " N = ", iN

            for iParticle in tupleParticles : # Loop over Particle types
                #print " Particle = ", iParticle
            
                # Create the ascii files needed for the energy resolution,
                # sampling fraction, and e/pi ratio plots, as a function
                # of the beam energy. These files, differently than those
                # for the shower shapes (longitudinal and transverse),
                # should not be created for each energy beam value. That's
                # why we have to create them here, rather than in the
                # function  extractInfo  , as instead for the shower shapes.
                # For the e/pi ratio, we consider only the particles  e-
                # and  pi+ .
                fileResolution_a = 0
                fileSampling_a   = 0
                fileRatio_a      = 0
                label_a = g4version1 + "-" + \
                          iPL + "-" + iCalo + "-" + iParticle + "-" + iN
                fileResolution_a = open( "energy_resolutions.txt-" + label_a , "w" )
                createdFiles.write( "energy_resolutions.txt-" + label_a + "\n" )
                fileSampling_a = open( "sampling_fractions.txt-" + label_a , "w" ) 
                createdFiles.write( "sampling_fractions.txt-" + label_a + "\n" )
                if ( iParticle == "e-"  or  iParticle == "pi+" ) :
                    fileRatio_a = open( "ratio_e_pi.txt-" + label_a , "w" )
                    createdFiles.write( "ratio_e_pi.txt-" + label_a + "\n" )
                fileResolution_b = 0
                fileSampling_b   = 0
                fileRatio_b      = 0
                if ( g4version2 ) :
                    label_b = g4version2 + "-" + \
                              iPL + "-" + iCalo + "-" + iParticle + "-" + iN
                    fileResolution_b = open( "energy_resolutions.txt-" + label_b , "w" )
                    createdFiles.write( "energy_resolutions.txt-" + label_b + "\n" )
                    fileSampling_b = open( "sampling_fractions.txt-" + label_b , "w" )
                    createdFiles.write( "sampling_fractions.txt-" + label_b + "\n" )
                    if ( iParticle == "e-"  or  iParticle == "pi+" ) :
                        fileRatio_b = open( "ratio_e_pi.txt-" + label_b , "w" )
                        createdFiles.write( "ratio_e_pi.txt-" + label_b + "\n" )

                for iE in tupleEnergies :     # Loop over Beam Energies
                    #print " Energy = ", iE
                    
                    foundFirst  = 0       # Found file for the first G4 version
                    foundSecond = 0       # Found file for the second G4 version
                    secondElement = 0     # Is it the second G4 version?
                    
                    for iG4 in tupleG4Versions : # Loop over Geant4 Versions
                        #print " G4 version = ", iG4

                        shortLABEL = iPL + "-" + iCalo + "-" \
                                   + iParticle + "-" + iE + "-" + iN
                        LABEL = iG4 + "-" + shortLABEL
                        #print " LABEL = ", LABEL
                        fileName = "output.log-" + LABEL
                        # Look for the file in all the subdirectories of
                        # the current directory.
                        command = "find " + directory + "/." + " -name " + fileName
                        command = command + " -follow -print > thePath.log"
                        #print '  command = ', command
                        os.system( command )

                        if ( os.path.getsize( "thePath.log" ) > 0 ) :
                            #print LABEL, " : FOUND! "
                            foundFiles.write( fileName + "\n" )

                            # Read the full path of the input file from
                            # thePath.log and then open the input file.
                            # In the unlikely case that there are more
                            # files with the same label, consider only the
                            # first one and neglect the others.
                            thePathFile = open( "thePath.log", "r" )
                            fullNameInputFile = ""
                            for line in thePathFile :
                                #print " line : ", line
                                fullNameInputFile = line.strip()
                                break
                            inputFile = open( fullNameInputFile, "r" )
                            if ( inputFile ) :                                
                                # Extract the useful information
                                if ( g4version2   and   secondElement ) :
                                    foundSecond = 1
                                    extractInfo( inputFile , LABEL , \
                                                 fileResolution_b , \
                                                 fileSampling_b , \
                                                 fileRatio_b )
                                else :
                                    foundFirst = 1
                                    extractInfo( inputFile , LABEL , \
                                                 fileResolution_a , \
                                                 fileSampling_a , \
                                                 fileRatio_a )
                            else :
                                print " inputFile = ", inputFile, " NOT FOUND!"

                            # Call the appropriate kumacs
                            if ( ( secondElement  or  len( tupleG4Versions ) == 1 ) and \
                                 ( foundFirst  or  foundSecond ) ) :
                                makeShowerPlots( shortLABEL, g4version1, g4version2 )
                                pass
                        else :
                            missingFiles.write( fileName + "\n" )
                            print LABEL, " : NOT found!"

                        secondElement = 1;

                # Close the ascii files needed for the energy resolution,
                # sampling fraction, and e/pi ratio plots (as a function
                # of the beam energy).
                fileResolution_a.close()
                fileSampling_a.close()
                if ( fileRatio_a ) :
                    fileRatio_a.close()
                if ( fileResolution_b ) :
                    fileResolution_b.close()
                if ( fileSampling_b ) :
                    fileSampling_b.close()
                if ( fileRatio_b ) :
                    fileRatio_b.close()

                aLABEL = iPL + "-" + iCalo + "-" + iParticle + "-" + iN
                makeResolutionPlots( aLABEL, g4version1, g4version2 )

            aLABEL = iPL + "-" + iCalo
            makeRatioPlots( aLABEL, iN, g4version1, g4version2 )
            
# Close the files.
foundFiles.close()
missingFiles.close()
createdFiles.close()
nonRetrievedFiles.close()

print '  ========== END observables.py ========== '


