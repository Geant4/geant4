#!/usr/bin/python

#-------------------------------------------------------------------
# Last update: 14-Jun-2006
#
# This script has no input arguments:
#
#         $  python createListPlots.py
#
# but assume that the text file, listPS.txt , which is created
# when the script  unpack.py  is run, is present in the directory
# where you run it.
# This script read the input file, and produced in output 3
# shell scripts:
#   1) plotAll.sh : which plots all the .ps files.
#   2) plotEnergies.sh : which plots only the visible energy
#                        (observable "1") and the total deposited
#                        energy (observable "2")
#   3) plotSelected.sh : which plots only the plots which match
#                        all the words (i.e. "and", not "or")
#                        specified in selectedCase
#                        (see ***LOOKHERE***).
#
# NB) After you run this script, you can use the generated shell
#     scripts as usual, e.g. :  $ ./plotAll.sh
#
#-------------------------------------------------------------------

import os
import sys
import string
import math

#***LOOKHERE***
selectedCases = (  "LHEP",
                   "WLAr",
                   "1GeV",
                   )


if ( len( sys.argv ) != 1 ) :
    print " Usage:  createListPlots.py "
else :

    outFileAll      = open( "plotAll.sh", "w" )
    outFileEnergies = open( "plotEnergies.sh", "w" )
    outFileSelected = open( "plotSelected.sh", "w" )

    theFile = open( "listPS.txt", "r" )
    countAll = 0
    countEnergies = 0
    countSelected = 0
    if ( theFile ) :
        for line in theFile :
            #print line
            countAll += 1
            outFileAll.write( "echo \" " + str( countAll ) + ")  " +
                              line.strip() + " \" ; gv " + line )
            if ( line.find( "-1-" ) > -1  or
                 line.find( "-2-" ) > -1 ) :
                countEnergies += 1
                outFileEnergies.write( "echo \" " + str( countEnergies ) + ")  " +
                                       line.strip() + " \" ; gv " + line )
            found = 1
            for word in selectedCases :
                if ( not line.find( word ) > -1 ) :
                    found = 0
            if ( found ) :
                countSelected += 1
                outFileSelected.write( "echo \" " + str( countSelected ) + ")  " +
                                       line.strip() + " \" ; gv " + line )
        print " countAll      = ", countAll
        print " countEnergies = ", countEnergies
        print " countSelected = ", countSelected
    else :
        print " No listPS.txt found! "
                
    theFile.close()
    outFileAll.close()
    outFileEnergies.close()
    outFileSelected.close()

    os.system( "chmod u+x plotAll.sh" )
    os.system( "chmod u+x plotEnergies.sh" )
    os.system( "chmod u+x plotSelected.sh" )
    

#-------------------------------------------------------------------

