#!/usr/bin/python

#--------------------------------------------------------------------------
# Last update: 18-Aug-2008
#
# This script requires at least 2 input arguments, and should be
# run as:
#
#         $  python addInfoInTable.py table file1 [file2] ... [fileN]
#
# where the  "table"  is a text (ascii) file with a table, where
# at least the first column, with the names of the observables,
# should be present, and the second argument, "file1" should be
# a log-files produced by running one of the  DetectorSlice
# simulations. The remain arguments,  file2 ... fileN  are optionals,
# and they are also meant to be log-files of  DetectorSlice  simulations.
# The table and files must be in the same directory as the script
# (use symbolic link if they are located in other directories).
#
# This script produces in output a new table, with the same name
# as the one given in input (whereas the original one is backup
# with name "old.0." in front of the original name; all other
# eventual intermediate tables are saved with names "old.1.",
# "old.2.", etc., in front of the original name: this is useful
# for debugging), in which new columns have been added with the
# corresponding information extract from the input files:
#    file1  first,  file2  second,  etc...
# using the script  printInfoLogFile.py  (which should be
# located in the same directory as addInfoInTable.py ).
# (For debugging, the output of  printInfoLogFile.py , when run
#  with each of the input files, can be found in:
#  .printInfoLogfile.out-file1 , .printInfoLogfile.out-file2, etc. )
#
# Notice that the rows corresponding to observables which are
# not recognized are reproduced in the new table as they are.
# If one of the files specified in the argument is not a log-file
# obtained by running one of the StatAccepTest simulations, nothing
# will happen, i.e. the table will not be updated.
#
# Look ***LOOKHERE*** to map the short names of the observables,
# used in the first column of the table, with the actual names
# of the observables as they appear in "file".
# However, if you have to build a new table, these are the names
# of the observables that can be used directly as the first column
# of the table, without requiring any change in the script below
# (of course, you can use a subset of them):
#
#        CPU
#        Tracker_E        
#        EM_Evis
#        EM_Etot
#        HAD_Evis
#        HAD_Etot
#        Muon_E
#        EM_res
#        HAD_res
#        EM_e_Evis
#        EM_mu_Evis
#        EM_pi_Evis
#        EM_k_Evis
#        EM_p_Evis
#        EM_n_Evis
#        HAD_e_Evis
#        HAD_mu_Evis
#        HAD_pi_Evis
#        HAD_k_Evis
#        HAD_p_Evis
#        HAD_n_Evis
#
# NB) The scripts energyScan.py , addInfoInTable.py , printInfoLogfile.py
#     have been obtained by modifying the Python scripts with the same
#     names in  StatAccepTest/dirEnergyScan/
#
#--------------------------------------------------------------------------

import os
import sys
import string
import math


#===============================================
#================= FUNCTIONS =================== 
#===============================================


def updateLine( infoList, tableLine ) :
    # This function received in input two arguments:
    # the first is a list whose elements are the lines of the
    # file with the information we want to add;
    # the second one is a line of the table which we want to
    # update with the information provided in the file.
    # This function returns the new, updated line.

    #***LOOKHERE*** : here is the map between the short names of the
    #                 observables (key of the dictionary), and the
    #                 complete name as it appears in "infoFile"
    #                 (value of the dictionary).
    mapNameObservables = {
        #
        'CPU':'cpuTime',
        'Tracker_E':'tracker_energy',
        'EM_Evis':'em_visEnergy',
        'EM_Etot':'em_totEnergy',
        'HAD_Evis':'had_visEnergy',
        'HAD_Etot':'had_totEnergy',
        'Muon_E':'muon_detector_energy',
        'EM_res':'em_resolution',
        'HAD_res':'had_resolution',
        'EM_e_Evis':'em_electronEvis',
        'EM_mu_Evis':'em_muonEvis',
        'EM_pi_Evis':'em_pionEvis',
        'EM_k_Evis':'em_kaonEvis',
        'EM_p_Evis':'em_protonEvis',
        'EM_n_Evis':'em_nucleusEvis',
        'HAD_e_Evis':'had_electronEvis',
        'HAD_mu_Evis':'had_muonEvis',
        'HAD_pi_Evis':'had_pionEvis',
        'HAD_k_Evis':'had_kaonEvis',
        'HAD_p_Evis':'had_protonEvis',
        'HAD_n_Evis':'had_nucleusEvis',
        #
        }

    #print ' tableLine=', tableLine

    newTableLine = tableLine

    #print ' tableLine.split()=', tableLine.split()
    if ( len( tableLine.split() ) > 0 ) :
        shortName = tableLine.split()[0]
        #print ' shortName=', shortName
        if ( mapNameObservables.has_key( shortName ) ) :
            longName = mapNameObservables[ shortName ]
            #print ' longName=', longName
            for infoLine in infoList :
                #print 'infoLine=', infoLine
                if ( infoLine.find( longName ) > -1 ) :
                    value = float( infoLine.split( "=" )[1] )
                    #print ' value=', value
                    if ( longName.find( "cpuTime" ) > -1 ) :
                         valueStringFormatted = "%.0f" % value
                         newTableLine += "   " + valueStringFormatted
                    else :
                         newTableLine += "   " + str( value )   
                         if ( longName.find( "resolution" ) > -1 ) :
                             newTableLine += "%"
                    break

    #print ' newTableLine=', newTableLine
    
    return newTableLine


#===============================================
#==================== MAIN ===================== 
#===============================================

# Check whether the number of arguments are correct or not.
if ( len( sys.argv ) < 3 ) :
    print " Usage:  addInfoInTable.py table file1 [file2] ... [fileN]"
else :

    # Loop over the input files: file1 , file2 ... fileN
    for iFile in xrange( len( sys.argv ) - 2 ) :
        print ' Considering file : ', sys.argv[ iFile+2 ]

        # Use the script  printInfoLogFile.py .
        os.system( "python printInfoLogfile.py " + sys.argv[ iFile+2 ] + \
                   " > .printInfoLogfile.out-" + sys.argv[ iFile+2 ] + " 2>&1 " )
    
        # The file  .printInfoLogfile.out-...  contains the
        # information we need in order to update the table.
        fileInformation = open( ".printInfoLogfile.out-" + sys.argv[ iFile+2 ], "r" )

        # Put the content of the file inside a list, in order to
        # be able to scan it more times.
        listInformation = []
        for infoLine in fileInformation :
            listInformation.append( infoLine.rstrip() )
            #print 'infoLine=', infoLine,

        print '  -> number of useful information = ', len( listInformation )
        #for lineInformation in listInformation :
        #    print 'lineInformation=', lineInformation
        
        fileInformation.close()

        if ( len( listInformation ) > 0 ) :
                
            table = open( sys.argv[1], "r" )
            if ( table ) :
                    
                # Create a new table.
                newtable = open( "new." + sys.argv[1], "w" )
                if ( newtable ) :
                    for line in table :
                        #print line,
                        newline = updateLine( listInformation, line.rstrip() )
                        newtable.write( newline + "\n" )

                newtable.close()
            table.close()

            # Rename the new table as the original one, but store 
            # also the previous table (for debugging).
            os.system( "mv " + sys.argv[1] + " old." + str( iFile ) + "." +
                       sys.argv[1] + 
                       " ; mv new." + sys.argv[1] + " " + sys.argv[1] )

#-------------------------------------------------------------------
