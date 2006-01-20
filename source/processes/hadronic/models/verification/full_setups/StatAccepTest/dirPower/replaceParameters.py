#!/usr/bin/python

#-------------------------------------------------------------------
# Last update: 20-Jan-2006
#
# This script received in input two parameters values: 
#
#         $  python  replaceParameters.py  alpha N
# e.g:
#         $  python  replaceParameters.py  0.5  3000
#
# which are interpreted as, respectively, the following two
# static constant parameters of  power.cpp :
#
#         1)  DistributionGenerator::alpha
#         2)  PowerCalculator::sampleSize 
#
# If the two specified input parameters have acceptable values,
# then the script modifies the program  power.cpp  by assigning
# those values to the corresponding static constants (the original
# program is saved as  save.power.cpp).
#
# In practice, this script is useful if you want to submit a bunch
# of jobs for power.cpp corresponding to different values of the
# the two above parameters.
# For instance, the shell script  jobs.sh  uses this Python script.
#
#-------------------------------------------------------------------

import sys
import os

###print " len( sys.argv ) = ", len( sys.argv )

# Check whether the number of arguments are correct or not.
isOK = 1
if ( len( sys.argv ) != 3 ) :
    print " Usage:  replaceParameters alpha N "
    isOK = 0
else:
    # Check whether the values of the arguments is acceptable or not.
    if ( float( sys.argv[1] ) < 0.0  or  float( sys.argv[1] ) > 1.0 ) :
        print " The parameter alpha must be between 0 and 1 : alpha= ", sys.argv[1]
        isOK = 0
    if ( int( sys.argv[2] ) <= 0  or  int( sys.argv[2] ) > 100000 ) :
        print " The parameter N must be non-negative and less than 100,000 : N= ", \
              sys.argv[2]
        isOK = 0

if ( isOK ) :
    # Open the file power.cpp , and create a new one.
    inFile = open( "power.cpp", 'r' )
    outFile = open( "temp", 'w' )

    # Scan the file inFile, line by line, and copy them to outFile,
    # but in the case the string "MAIN-PARAMETER" appears, then
    # assign the corresponding parameter before copying it to outFile.
    for line in inFile :
        if ( line.find( "MAIN-PARAMETER" ) >= 0 ) :
            newline = ""
            if ( line.find( "alpha" ) >= 0 ) :
                newline = "const double DistributionGenerator::alpha = " + sys.argv[1]  
            elif ( line.find( "sampleSize" ) >= 0 ) :
                newline = "const int PowerCalculator::sampleSize       = " + sys.argv[2]
            newline += ";   //***MAIN-PARAMETER*** \n" 
            outFile.write( newline )
        else :
            outFile.write( line )
    
    # Close the files.
    inFile.close()
    outFile.close()

    # Rename the files:  power.cpp  -> save.power.cpp
    #                    temp -> power.cpp
    os.rename( "power.cpp", "save.power.cpp" )
    os.rename( "temp", "power.cpp" )

#-------------------------------------------------------------------

