#!/usr/bin/python

#-----------------------------------------------------------------
# Last update: 05-Aug-2008
# 
# This Python script has one input argument, which is a Geant4
# command file for the simplified calorimeter, e.g.
#
#        python writeFlukaCaloGeometry.py atlasHEC.g4
#
# and it outputs 2 files:
#
#    1)  mysmallinclude.inc : a Fortran include file to be
#                             copied in the subdirectory
#                             dirUserRoutines/ .
#
#    2)  caloInput-xxx.inp : the Fluka input file, corresponding
#                            to the simplified calorimeter
#                            specified in the Geant4 command file.
#                            The string "xxx" identifies uniquely
#                            the configuration (beam particle,
#                            beam energy, calorimeter).
#
# Notice that:
#
#   -  To avoid the complications of writing with a format,
#      THE FILE  caloInput-xxx.inp  NEEDS TO BE EDITED IN ORDER
#      TO PROPERLY FIX THE EXPECTED POSITIONS OF THE VARIOUS
#      FIELDS. This is quite straightforward, because templates
#      for the expected fields are present in the file.
#      But if you do not do it, it will not work!
#
#   -  Magnetic field is ignored (in Fluka);
#
#   -  The kinetic energy thresholds are corresponding to the
#      default production range cut of Geant4, i.e. 0.7 mm.
#      See "***LOOKHERE***THRESHOLDS" to change them.
#
#   -  Look for the string "***LOOKHERE***BIRKS" to set the
#      Birks coefficients, in the case you want to use Birks
#      quenching (you need also to edit the file mgdraw.f ...).
#
#   -  Look for the string "***LOOKHERE***" for all the places
#      of this script where parameters are hardwired...
#
#-----------------------------------------------------------------

import os
import sys
import string
import math

#===============================================
#================= GLOBAL PARAMETERS ===========
#===============================================

#***LOOKHERE*** : if the following variable is "0"
#                 then the new style (introduced in
#                 September 2007) for the material
#                 definitions is used;
#                 else, the old style is used.
isOldStyleMaterialDefinitions = 0 


#===============================================
#================= FUNCTIONS =================== 
#===============================================

# ====== Some functions for string manipulation ======

def funExtract( theString, theWord ) :
    # Given in input the string "theString", this function
    # search for "theWord", and if finds it, it returns the
    # substring of "theString" which is after "theWord".
    # For example:  funExtract( "value 7.13", "value" )
    # returns  " 7.13" .
    result = ""
    if ( theString.find( theWord ) >= 0 ) :
        substring = theString.split( theWord )[1]
        #print ' substring=', substring
        pos = 0
        while ( pos < len( substring )  and  not substring[pos].isalnum() ) :
            pos += 1
        if ( pos < len( substring ) ) :
            result = substring[pos:]
        #print ' pos=', pos, '  result=', result
    return result


def writeFunction( beamParticleInput, beamEnergyInput,
                   absorberMaterialInput, activeMaterialInput,
                   isCalHomogeneousInput, isUnitInLambdaInput,
                   absorberTotalLengthInput, calorimeterRadiusInput,
                   activeLayerNumberInput, readoutLayerNumberInput,
                   activeLayerSizeInput, isRadiusUnitInLambdaInput,
                   radiusBinSizeInput, radiusBinNumberInput,
                   numberEventsInput ) :
    # This is the main function that writes the two output files.
    # Checks are made to ensure that the input arguments have
    # their expected values.
    # It is essential to properly convert  string -> numbers ,
    # when you need to manipulate numbers, and vice versa,
    # convert  numbers -> string  when you need to output (i.e.
    # write into a file) those numbers.

    # ---------------- Particle type ---------------
    particleType = ""
    dictParticle = { 'mu-':'MUON-' , 'mu+':'MUON+' ,
                     'e-':'ELECTRON' , 'e+':'POSITRON' , 'gamma':'PHOTON' ,
                     'pi+':'PION+'   , 'pi-':'PION-'  ,
                     'kaon+':'KAON+'  , 'kaon-':'KAON-' , 'kaon0L':'KAONLONG' ,
                     'neutron':'NEUTRON' , 'proton':'PROTON' ,
                     'anti_neutron':'ANEUTRON' , 'anti_proton':'APROTON' }
    if ( dictParticle.has_key( beamParticleInput ) ) :
        particleType = dictParticle[ beamParticleInput ]
    else :
        print '  ***ERROR*** in writeFunction.py : WRONG beamParticleInput = ', \
              beamParticleInput
        sys.exit(0)        
    print '  particleType = ', particleType
                
    # ---------------- Beam energy -----------------
    energyValue = ""
    numericEnergyValue = 0.0
    isNumericPart = 1
    for character in beamEnergyInput :
        if ( isNumericPart ) :
             if ( character.isdigit() ) :
                energyValue = energyValue + character
                numericEnergyValue = float( energyValue )
             elif ( character.isalpha() ) :
                numericEnergyValue = float( energyValue )
                energyValue = energyValue + " " + character
                isNumericPart = 0
        else :
            if ( character.isalpha() ) :
                energyValue = energyValue + character
            elif ( character.isdigit() ) :
                print '  ***ERROR*** in writeFunction.py : WRONG beamEnergyInput = ', \
                      beamEnergyInput
                sys.exit(0)
    # We need the beam kinetic energy in GeV : by default, if unit
    # is not specified, the assumed unit in Geant4 is MeV; we also
    # assume it if you do not understand the unit.
    if ( energyValue.find( " GeV" ) >= 0 ) :
        pass
    elif ( energyValue.find( " MeV" ) >= 0 ) :
        numericEnergyValue /= 1000.0
    elif ( energyValue.find( " keV" ) >= 0 ) :
        numericEnergyValue /= 1000000.0
    elif ( energyValue.find( " TeV" ) >= 0 ) :
        numericEnergyValue *= 1000.0
    elif ( energyValue.find( " eV" ) >= 0 ) :
        numericEnergyValue /= 1000000000.0
    else :
        print '-> NOT recognized beam energy unit, MeV is ASSUMED '
        numericEnergyValue /= 1000.0
    print '  energyValue = ', energyValue, \
          '  numericEnergyValue = ', numericEnergyValue 

    # ---------------- Calorimeter properties ------------
    absorberMaterial = ""                                 #***LOOKHERE***ABSORBER
    dictMaterial = { 'Fe':'IRON' ,    'Iron':'IRON',
                     'Cu':'COPPER',   'Copper':'COPPER',
                     'W':'TUNGSTEN',  'Tungsten':'TUNGSTEN', 
                     'Pb':'LEAD',     'Lead':'LEAD',
                     'U':'URANIUM',   'Uranium':'URANIUM',
                     'PbWO4':'PBWO4' }
                     
    if ( dictMaterial.has_key( absorberMaterialInput ) ) :
        absorberMaterial = dictMaterial[ absorberMaterialInput ]
    else :
        for fullName in dictMaterial.values() :
            if ( fullName == absorberMaterialInput  or
                 fullName.lower() == absorberMaterialInput  or
                 fullName.capitalize() == absorberMaterialInput  or
                 fullName.lower().swapcase() == absorberMaterialInput ) :
                absorberMaterial = fullName
            if ( len( absorberMaterial ) == 0 ) :
                print '  ***ERROR*** in writeFunction.py : WRONG absorberMaterialInput = ', absorberMaterialInput
                sys.exit(0)        
    print '  absorberMaterial = ', absorberMaterial

    activeMaterial = ""                                   #***LOOKHERE***ACTIVE MAT.
    dictMaterial = { 'Scintillator':'POLYSTYR' ,    'Sci':'POLYSTYR',
                     'LiquidArgon':'LIQUIDAR',   'LAr':'LIQUIDAR',
                     'PbWO4':'PBWO4' ,
                     'Silicon':'SILICON', 'Si':'SILICON' }
                     
    if ( dictMaterial.has_key( activeMaterialInput ) ) :
        activeMaterial = dictMaterial[ activeMaterialInput ]
    else :
        for fullName in dictMaterial.values() :
            if ( fullName == activeMaterialInput  or
                 fullName.lower() == activeMaterialInput  or
                 fullName.capitalize() == activeMaterialInput  or
                 fullName.lower().swapcase() == activeMaterialInput ) :
                activeMaterial = fullName
            if ( len( activeMaterial ) == 0 ) :
                print '  ***ERROR*** in writeFunction.py : WRONG activeMaterialInput = ', activeMaterialInput
                sys.exit(0)        
    print '  activeMaterial = ', activeMaterial

    # Interaction lengths of different materials, in [cm]
    dictLambda = { '':0.0,                                #***LOOKHERE***LAMBDA VALUES
                   'IRON':16.760,
                   'COPPER':15.056,
                   'TUNGSTEN':9.5855,
                   'LEAD':17.092, 
                   'URANIUM':10.501,
                   'PBWO4':22.4 }

    absorberTotalLengthCM = 0.0   # In centimeters.
    calorimeterRadiusCM = 0.0     #  "     "
    if ( isUnitInLambdaInput.find( "1" ) >= 0 ) :
        absorberTotalLengthCM = float( absorberTotalLengthInput ) *     \
                                dictLambda[ absorberMaterial ]
        calorimeterRadiusCM = float( calorimeterRadiusInput ) *         \
                              dictLambda[ absorberMaterial ]
    else :
        absorberTotalLengthCM = float( absorberTotalLengthInput ) / 10.0
        calorimeterRadiusCM = float( calorimeterRadiusInput ) / 10.0
    print '  absorberTotalLengthCM = ', absorberTotalLengthCM, ' [cm]'
    print '  calorimeterRadiusCM   = ', calorimeterRadiusCM, ' [cm]'

    activeLayerNumber = int( activeLayerNumberInput )
    print '  activeLayerNumber = ', activeLayerNumber
    
    readoutLayerNumber = int( readoutLayerNumberInput )
    print '  readoutLayerNumber = ', readoutLayerNumber

    numberOfLayersPerReadoutLayer = activeLayerNumber / readoutLayerNumber
    if ( math.fmod( activeLayerNumber, readoutLayerNumber ) != 0 ) :
        print '  numberOfLayersPerReadoutLayer : NOT AN INTEGER!'
    else :
        print '  numberOfLayersPerReadoutLayer = ', numberOfLayersPerReadoutLayer
    
    activeLayerSizeCM = float( activeLayerSizeInput ) / 10.0  # In centimeters.
    print '  activeLayerSizeCM = ', activeLayerSizeCM, ' [cm]'

    totalLengthCalorimeterCM = absorberTotalLengthCM
    thicknessAbsorberLayerCM = absorberTotalLengthCM / float( activeLayerNumber )
    if ( isCalHomogeneousInput.find( "0" ) >= 0 ) :
        totalLengthCalorimeterCM += activeLayerNumber * activeLayerSizeCM
    else :
        thicknessAbsorberLayerCM = ( absorberTotalLengthCM - \
                                     activeLayerNumber * activeLayerSizeCM ) \
                                     / float( activeLayerNumber )
    print '  totalLengthCalorimeterCM = ', totalLengthCalorimeterCM, ' [cm]'
    print '  thicknessAbsorberLayerCM = ', thicknessAbsorberLayerCM, ' [cm]'

    radiusBinSizeCM = 0.0   # In centimeters.
    if ( isRadiusUnitInLambdaInput.find( "1" ) >= 0 ) :
        radiusBinSizeCM = float( radiusBinSizeInput ) * dictLambda[ absorberMaterial ]
    else :
        radiusBinSizeCM = float( radiusBinSizeInput ) / 10.0
    print '  radiusBinSizeCM = ', radiusBinSizeCM, ' [cm]'

    radiusBinNumber = int( radiusBinNumberInput )
    print '  radiusBinNumber = ', radiusBinNumber

    # ----------------- Write the Fortran include file  ----------

    fileInclude = open( "mysmallinclude.inc", "w" )
    
    fileInclude.write( "      INTEGER maxNumberOfEvents\n" )
    fileInclude.write( "      INTEGER numberOfReadoutLayers\n" )
    fileInclude.write( "      INTEGER numberOfLayersPerReadoutLayer\n" )
    fileInclude.write( "      INTEGER numberOfRadiusBins\n" )
    fileInclude.write( "      REAL radiusBin\n" )

    #***LOOKHERE*** maxNumberOfEvents value 
    fileInclude.write( "      PARAMETER ( maxNumberOfEvents = 10000 )\n" )

    fileInclude.write( "      PARAMETER ( numberOfReadoutLayers = " + \
                       readoutLayerNumberInput + " )\n" )
    fileInclude.write( "      PARAMETER ( numberOfLayersPerReadoutLayer = " + \
                       str( numberOfLayersPerReadoutLayer ) + " )\n" )
    fileInclude.write( "      PARAMETER ( numberOfRadiusBins = " + \
                       radiusBinNumberInput + " )\n" )
    fileInclude.write( "      PARAMETER ( radiusBin = " + \
                       str( radiusBinSizeCM ) + " )\n" )

    fileInclude.close()

    # ----------------- Write Fluka input file -------------------

    fileFlukaInput = "caloInput-" + particleType + "-" + \
                     str( int( numericEnergyValue ) ) + "GeV-" + \
                     absorberMaterial + "-" + activeMaterial + "-" + \
                     numberEventsInput + "evts"
    
    theFile = open( fileFlukaInput + ".inp" , "w" )
    
    theFile.write( "TITLE\n" )
    theFile.write( "  Simplified calorimeter\n" )
    theFile.write( "**COD+NUM+---1st----+++2nd++++---3rd----+++4th++++---5th----+++6th++++\n" )

    if ( isOldStyleMaterialDefinitions ) :
        theFile.write( "GLOBAL                                         4.0\n" )


    theFile.write( "PHYSICS        100.0     100.0     100.0     100.0     100.0     100.0PEATHRES\n" )
    
    theFile.write( "DEFAULTS                                                              CALORIME\n" )
    
    # A negative value means kinetic energy for the beam, instead of momentum.
    theFile.write( "BEAM         -" + str( numericEnergyValue ) + \
                   "                                                   " + \
                   particleType + "\n" )

    #***LOOKHERE***BEAM POSITION
    theFile.write( "BEAMPOS          0.0       0.0    -200.0\n" )
    
    theFile.write( "GEOBEGIN                  0.01                                        COMBINAT\n" )
    theFile.write( "                    simple calorimeter\n" )
    theFile.write( "**COD+NUM+---1st----+++2nd++++---3rd----+++4th++++---5th----+++6th++++\n" )

    # Start of the Body data.
    
    #***LOOKHERE***EXPERIMENTAL HALL
    theFile.write( "  RPP    1   -2000.0   +2000.0   -2000.0   +2000.0   -2000.0   +2000.0\n" )
    theFile.write( "  RPP    2   -1000.0   +1000.0   -1000.0   +1000.0   -1000.0   +1000.0\n" )

    # We use an "Infinite Circular Cylinder parallel to z coordinate axis"
    # (ZCC), centered at (x=0.0, y=0.0) with a specified radius, and then
    # we limit (or slice) it by using a set of "Infinite half-spaces 
    # delimited by z coordinate plane" (XYP), in correspondence of the
    # starting of each layer (passive or active).
    theFile.write( "  ZCC    3      0.00       0.0 " + \
                   str( calorimeterRadiusCM ) + "\n" )

    zPosition = - totalLengthCalorimeterCM / 2.0
    for i in xrange( 2*activeLayerNumber + 1 ) :
        theFile.write( "  XYP " + str( i+4 ) + "   " + str( zPosition ) + "\n" )
        if ( math.fmod( i, 2 ) == 0 ) :
            zPosition += thicknessAbsorberLayerCM 
        else :
            zPosition += activeLayerSizeCM
        if ( math.fabs( zPosition ) < 1.0e-6 ) :
            zPosition = 0.0

    theFile.write( "  END\n" )

    # Start of the Region data
    theFile.write( "**REG     cc1st--cc2nd--cc3rd--cc4th--cc5th--cc6th--cc7th--cc8th--cc9th--\n" )

    # Black hole
    theFile.write( "    1          +1     -2\n" )

    # Gap (of air) in between black hole and the calorimeter
    # the first zone is the entire experimental hall after removing the
    # infinite cylinder; the second zone is the cylinder of air in 
    # between the experimental hall and the calorimeter, at z < 0;
    # the third zone is the cylinder of air in betwen the calorimeter
    # and the experimental hall, at z > 0. 
    numberOfRegions = 4 + 2*activeLayerNumber
    theFile.write( "    2     OR   +2     -3OR   +3     +4     +2OR   +3   -" + \
                   str( numberOfRegions ) + "     +2\n" )

    # Now describe the calorimeter, in a flat way.
    # Notice that:  odd  regions  are passive layers (absorber);
    #               even regions  are active  layers.   
    for i in xrange( 2*activeLayerNumber ) :
        theFile.write( "  " + str( i+3 ) + "          +3   -" + \
                       str( i+4 ) + "     +" + str( i+5 ) + "\n" )
    
    theFile.write( "  END\n" )
    theFile.write( "GEOEND\n" )
    theFile.write( "****--- Debug the geometry\n" )
    theFile.write( "***GEOEND        1000.0    1000.0    1000.0   -1000.0   -1000.0   -1000.0DEBUG\n" )
    theFile.write( "***GEOEND         200.0     200.0     200.0                              &\n" )
    theFile.write( "*-CODEWD  ---1st----+++2nd++++---3rd----+++4th++++---5th----+++6th++++SDUM----\n" )

    #***LOOKHERE***MATERIAL PROPERTIES
    if ( isOldStyleMaterialDefinitions ) :
        theFile.write( "MATERIAL        26.0    55.850 7.870e+00       3.0                    IRON\n" )
        theFile.write( "MATERIAL        29.0    63.540 8.960e+00       4.0                    COPPER\n" )
        theFile.write( "MATERIAL        74.0   183.850 1.930e+01       5.0                    TUNGSTEN\n" )
        theFile.write( "MATERIAL        82.0   207.190 1.135e+01       6.0                    LEAD\n" )
        theFile.write( "MATERIAL        92.0   238.030 1.895e+01       7.0                    URANIUM\n" )
        theFile.write( "MATERIAL        18.0    39.950 1.400e+00       8.0                    LIQUIDAR\n" )
        theFile.write( "MATERIAL                       1.290e-03       9.0                    AIR\n" )
        theFile.write( "MATERIAL         7.0    14.010 9.990e-01      10.0                    NITROGEN\n" )
        theFile.write( "MATERIAL         8.0    16.000 9.990e-01      11.0                    OXYGEN\n" )
        theFile.write( "MATERIAL                       1.000e-05      12.0                    G4VACUUM\n" )
        theFile.write( "MATERIAL                       1.032e+00      13.0                    POLYSTYR\n" )
        theFile.write( "MATERIAL         6.0    12.010 9.990e-01      14.0                    CARBON\n" )
        theFile.write( "MATERIAL         1.0     1.010 9.990e-01      15.0                    HYDROGEN\n" )
        theFile.write( "MATERIAL                       8.280e+00      16.0                    PBWO4\n" )
        theFile.write( "MATERIAL        14.0    28.085 2.330e+00      17.0                    SILICON\n" )
        theFile.write( "LOW-MAT          8.0      18.0      -2.0      87.0       0.0      0.0 ARGON\n" )
    else :
        theFile.write( "MATERIAL        26.0    55.850 7.870e+00       3.0                    MY_FE\n" )
        theFile.write( "MATERIAL        29.0    63.540 8.960e+00       4.0                    MY_CU\n" )
        theFile.write( "MATERIAL        74.0   183.850 1.930e+01       5.0                    MY_W\n" )
        theFile.write( "MATERIAL        82.0   207.190 1.135e+01       6.0                    MY_PB\n" )
        theFile.write( "MATERIAL        92.0   238.030 1.895e+01       7.0                    MY_U\n" )
        theFile.write( "MATERIAL        18.0    39.950 1.400e+00       8.0                    LIQUIDAR\n" )
        theFile.write( "MATERIAL                       1.290e-03       9.0                    AIR\n" )
        theFile.write( "MATERIAL         7.0    14.010 9.990e-01      10.0                    MY_N\n" )
        theFile.write( "MATERIAL         8.0    16.000 9.990e-01      11.0                    MY_O\n" )
        theFile.write( "MATERIAL                       1.000e-05      12.0                    G4VACUUM\n" )
        theFile.write( "MATERIAL                       1.032e+00      13.0                    POLYSTYR\n" )
        theFile.write( "MATERIAL         6.0    12.010 9.990e-01      14.0                    MY_C\n" )
        theFile.write( "MATERIAL         1.0     1.010 9.990e-01      15.0                    MY_H\n" )
        theFile.write( "MATERIAL                       8.280e+00      16.0                    PBWO4\n" )
        theFile.write( "MATERIAL        14.0    28.085 2.330e+00      17.0                    MY_SI\n" )
        theFile.write( "LOW-MAT          3.0                                                  IRON\n" )
        theFile.write( "LOW-MAT          4.0                                                  COPPER\n" )
        theFile.write( "LOW-MAT          5.0                                                  TUNGSTEN\n" )
        theFile.write( "LOW-MAT          6.0                                                  LEAD\n" )
        theFile.write( "LOW-MAT          7.0                                                  238-U\n" )
        theFile.write( "LOW-MAT          8.0      18.0      -2.0      87.0       0.0      0.0 ARGON\n" )
        theFile.write( "LOW-MAT         10.0                                                  NITROGEN\n" )
        theFile.write( "LOW-MAT         11.0                                                  OXYGEN\n" )
        theFile.write( "LOW-MAT         14.0                                                  CARBON\n" )
        theFile.write( "LOW-MAT         15.0                                                  HYDROGEN\n" )
        theFile.write( "LOW-MAT         17.0                                                  SILICON\n" )

    # Air compound
    theFile.write( "COMPOUND   -0.700000      10.0 -0.300000      11.0                    AIR\n" )

    # Geant4 Vacuum compound (i.e. air at lower pressure)
    theFile.write( "COMPOUND   -0.700000      10.0 -0.300000      11.0                    G4VACUUM\n" )

    # PbWO4 compound
    theFile.write( "COMPOUND   -0.455323       6.0 -0.404030       5.0 -0.140647     11.0 PBWO4\n" )

    # Scintillator compound
    theFile.write( "COMPOUND   -0.914956      14.0 -0.085044      15.0                    POLYSTYR\n" )
    
    theFile.write( "*-CODEWD  ---1st----+++2nd++++---3rd----+++4th++++---5th----+++6th++++SDUM----\n" )

    # External black hole
    theFile.write( "ASSIGNMAT        1.0       1.0\n" )

    # Gap: assign the Geant4 vacuum (i.e. air at lower pressure)
    theFile.write( "ASSIGNMAT       12.0       2.0\n" )

    absorberMaterialNum = 0
    cutElectronsString = 0.0
    cutPhotonsString = 0.0                   #***LOOKHERE***THRESHOLDS in GeV
    if ( absorberMaterial == "IRON" ) :
        absorberMaterialNum = 3
        cutElectronsString = "0.953E-3"
        cutPhotonsString = "17.1E-6"
    elif ( absorberMaterial == "COPPER" ) :
        absorberMaterialNum = 4
        cutElectronsString = "1.026E-3"
        cutPhotonsString = "20.5E-6"
    elif ( absorberMaterial == "TUNGSTEN" ) :
        absorberMaterialNum = 5
        cutElectronsString = "1.637E-3"
        cutPhotonsString = "97.4E-6"
    elif ( absorberMaterial == "LEAD" ) :
        absorberMaterialNum = 6
        cutElectronsString = "1.001E-3"
        cutPhotonsString = "94.7E-6"
    elif ( absorberMaterial == "URANIUM" ) :
        absorberMaterialNum = 7
        cutElectronsString = "1.521E-3"
        cutPhotonsString = "113.4E-6"
    elif ( absorberMaterial == "PBWO4" ) :
        absorberMaterialNum = 16
        cutElectronsString = "0.842E-3"
        cutPhotonsString = "74.7E-6"
    else :
        print '  NOT RECOGNIZED absorberMaterial = ', absorberMaterial

    activeMaterialNum = 0
    if ( activeMaterial == "LIQUIDAR" ) :
        activeMaterialNum = 8
    elif ( activeMaterial == "POLYSTYR" ) :
        activeMaterialNum = 13
    elif ( activeMaterial == "PBWO4" ) :
        activeMaterialNum = 16
    elif ( activeMaterial == "SILICON" ) :
        activeMaterialNum = 17
    else :
        print '  NOT RECOGNIZED activeMaterial = ', activeMaterial

    # Calorimeter : Odd  regions: passive layers (absorber); 
    #               Even regions: active layers.
    for i in xrange( numberOfRegions-4 ) :
        materialNum = 0
        if ( math.fmod( i, 2 ) == 0 ) :
            materialNum = absorberMaterialNum
        else :
            materialNum = activeMaterialNum
        theFile.write( "ASSIGNMAT       " + str( materialNum ) + \
                       ".0     " + str( i+3 ) + ".0\n" )
                   
    theFile.write( "*-CODEWD  ---1st----+++2nd++++---3rd----+++4th++++---5th----+++6th++++SDUM----\n" )

    # Physics

    #***LOOKHERE***PHYSICS OPTIONS

    # To switch off muon-nuclear process. It is on by default in CALORIME.
    ###theFile.write( "***MUPHOTON        -1.0\n" )

    # To switch on gamma-nuclear process. It is off by default in CALORIME.
    theFile.write( "PHOTONUC         1.0       0.0       0.0       3.0      99.0\n" )

    # To switch off the low-energy neutron transport. It is on by default in CALORIME.
    theFile.write( "***LOW-BIAS         1.0                                   999.0\n" )

    # Thresholds: electrons/positrons/photons production and tranportation cutoffs.
    theFile.write( "EMFCUT     -" + cutElectronsString + \
                   "   " + cutPhotonsString + \
                   "       1.0       3.0      99.0          PROD-CUT\n" )

    # Muons/hadrons tranportation cutoff.
    # Default for CALORIME is  (kinetic energy) 0.001*m/m_p GeV .
    # There is no transportation cut off in Geant4, so the best
    # thing to do is to set it to the minimimum possible in Fluka,
    # which seems to be 100 keV. 
    ###theFile.write( "PART-THR     -0.0001       1.0      62.0                 0.0\n" )

    # Muons/hadrons bremsstrahlung (pair production always allowed).
    # Default for CALORIME is 300 keV. Set here at 100 keV.
    ###theFile.write( "PAIRBREM                          0.0001       3.0      99.0\n" )

    #***endLOOKHERE*** (physics options)

    # Activate transportation of neutrinos
    theFile.write( "DISCARD         -5.0      -6.0     -27.0     -28.0     -43.0     -44.0\n" )

    theFile.write( "SCORE          201.0     208.0\n" )

    # Activate the run initialization :  USRINI .
    theFile.write( "USRICALL\n" )

    # Activate the "stepping action" : MGDRAW .
    theFile.write( "USERDUMP       200.0                           1.0\n" )

    #***LOOKHERE***BIRKS : select the Birks coefficients
    #
    # ATLAS & CMS scintillator tiles (rho = 1.032 [gr/cm^3] )
    ###theFile.write( "USERDUMP     1.29E-2   9.59E-6                                        UDQUENCH\n" )     
    #
    # Standard Birks scintillator coefficients (rho = 1.032 [gr/cm^3] )
    ###theFile.write( "USERDUMP     1.31E-2       0.0                                        UDQUENCH\n" )     
    #
    # Zeus scintillator SCSN38 (lower limit) (rho = 1.032 [gr/cm^3] )
    ###theFile.write( "USERDUMP      8.5E-3       0.0                                        UDQUENCH\n" )     
    #
    # Pilot B scintillator (same paper, upper limit) (rho = 1.032 [gr/cm^3] ) 
    ###theFile.write( "USERDUMP     1.59E-2       0.0                                        UDQUENCH\n" )     
    #
    # Liquid Argon : my estimation for ATLAS ( rho = 1.396 [gr/cm^3] )
    ###theFile.write( "USERDUMP       0.022       0.0                                        UDQUENCH\n" )     
    #
    # Liquid Argon : most used value
    ###theFile.write( "USERDUMP       0.005       0.0                                        UDQUENCH\n" )
    #
    #*********************************

    # Select random seed  
    theFile.write( "RANDOMIZ         1.0 7691919.0\n" )   #***LOOKHERE***RANDOM SEED

    # Select the number of events.
    theFile.write( "START         " + numberEventsInput + ".0            999999.0                 0.0\n" )

    # Activate the run end :  USROUT .
    theFile.write( "USROCALL\n" )
    
    theFile.write( "STOP\n" )
     
    theFile.close()

    return 0


#===============================================
#==================== MAIN ===================== 
#===============================================

print '========== START writeFlukaCaloGeometry.py =========='

###print ' len( sys.argv ) = ', len( sys.argv )

# Check whether the number of arguments are correct or not.
if ( len( sys.argv ) != 2 ) :
    print " Usage:  writeFlukaCaloGeometry file1 "
else :
    # Open the file given in input.
    file1 = sys.argv[1]

    print ' '
    print ' input file  : ', file1
    print ' '
 
    inFile1 = open( file1, 'r' )

    beamParticleString = ""
    beamEnergyString = ""
    absorberMaterialString = ""
    activeMaterialString = ""
    isCalHomogeneousString = ""
    isUnitInLambdaString = ""
    absorberTotalLengthString = ""
    calorimeterRadiusString = ""
    activeLayerNumberString = ""
    readoutLayerNumberString = ""
    activeLayerSizeString = ""
    isRadiusUnitInLambdaString = ""
    radiusBinSizeString = ""
    radiusBinNumberString = ""
    numberEventsString = ""

    for line in inFile1 :
        #print ' line=', line.rstrip()
        pos = 0
        while ( pos < len( line.rstrip() )  and  line[ pos ] != "#"  ) :
            pos += 1
        if ( pos > 0 ) :
            rline = line.rstrip()[0:pos]
            #print ' rline=', rline
            if ( rline.find( "/gun/particle" ) >= 0 ) :
                beamParticleString = funExtract( rline, "/gun/particle" )
            if ( rline.find( "/gun/energy" ) >= 0 ) :
                beamEnergyString = funExtract( rline, "/gun/energy" )
            if ( rline.find( "/mydet/absorberMaterial" ) >= 0 ) :
                absorberMaterialString = funExtract( rline, "/mydet/absorberMaterial" )
            if ( rline.find( "/mydet/activeMaterial" ) >= 0 ) :
                activeMaterialString = funExtract( rline, "/mydet/activeMaterial" )
            if ( rline.find( "/mydet/isCalHomogeneous" ) >= 0 ) :
                isCalHomogeneousString = funExtract( rline, "/mydet/isCalHomogeneous" )
            if ( rline.find( "/mydet/isUnitInLambda" ) >= 0 ) :
                isUnitInLambdaString = funExtract( rline, "/mydet/isUnitInLambda" )
            if ( rline.find( "/mydet/absorberTotalLength" ) >= 0 ) :
                absorberTotalLengthString = funExtract( rline,
                                                        "/mydet/absorberTotalLength" )
            if ( rline.find( "/mydet/calorimeterRadius" ) >= 0 ) :
                calorimeterRadiusString = funExtract( rline, "/mydet/calorimeterRadius" )
            if ( rline.find( "/mydet/activeLayerNumber" ) >= 0 ) :
                activeLayerNumberString = funExtract( rline, "/mydet/activeLayerNumber" )
            if ( rline.find( "/mydet/readoutLayerNumber" ) >= 0 ) :
                readoutLayerNumberString = funExtract( rline,
                                                       "/mydet/readoutLayerNumber" )
            if ( rline.find( "/mydet/activeLayerSize" ) >= 0 ) :
                activeLayerSizeString = funExtract( rline, "/mydet/activeLayerSize" )
            if ( rline.find( "/mydet/isRadiusUnitInLambda" ) >= 0 ) :
                isRadiusUnitInLambdaString = funExtract( rline,
                                                         "/mydet/isRadiusUnitInLambda" )
            if ( rline.find( "/mydet/radiusBinSize" ) >= 0 ) :
                radiusBinSizeString = funExtract( rline, "/mydet/radiusBinSize" )
            if ( rline.find( "/mydet/radiusBinNumber" ) >= 0 ) :
                radiusBinNumberString = funExtract( rline, "/mydet/radiusBinNumber" )
            if ( rline.find( "/run/beamOn" ) >= 0 ) :
                numberEventsString = funExtract( rline, "/run/beamOn" )

            if ( rline.find( "/mydet/setField" ) >= 0 ) :
                print ' WARNING : MAGNETIC FIELD IS IGNORED! ', rline

    if ( not readoutLayerNumberString ) :
        readoutLayerNumberString = activeLayerNumberString
        print ' WARNING : readoutLayerNumber NOT SPECIFIED! ASSUMED TO BE = activeLayerNumber'
    print ' beamParticleString          =', beamParticleString        
    print ' beamEnergyString            =', beamEnergyString        
    print ' absorberMaterialString      =', absorberMaterialString
    print ' activeMaterialString        =', activeMaterialString
    print ' isCalHomogeneousString      =', isCalHomogeneousString
    print ' isUnitInLambdaString        =', isUnitInLambdaString
    print ' absorberTotalLengthString   =', absorberTotalLengthString
    print ' calorimeterRadiusString     =', calorimeterRadiusString
    print ' activeLayerNumberString     =', activeLayerNumberString
    print ' readoutLayerNumberString    =', readoutLayerNumberString
    print ' activeLayerSizeString       =', activeLayerSizeString
    print ' isRadiusUnitInLambdaString  =', isRadiusUnitInLambdaString
    print ' radiusBinSizeString         =', radiusBinSizeString
    print ' radiusBinNumberString       =', radiusBinNumberString
    print ' numberEventsString          =', numberEventsString        
    print ' '
    
    writeFunction( beamParticleString, beamEnergyString,
                   absorberMaterialString, activeMaterialString,
                   isCalHomogeneousString, isUnitInLambdaString,
                   absorberTotalLengthString, calorimeterRadiusString,
                   activeLayerNumberString, readoutLayerNumberString,
                   activeLayerSizeString, isRadiusUnitInLambdaString,
                   radiusBinSizeString, radiusBinNumberString,
                   numberEventsString )

    # Close the files.
    inFile1.close()

print ' '
print '========== END writeFlukaCaloGeometry.py =========='
