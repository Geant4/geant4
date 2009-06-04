#!/usr/bin/python

#--------------------------------------------------------------------------
# Last update: 04-Jun-2009
#
# This script requires no input argument, and should be run as:
#
#         $  python runningSites.py
#
# The script examines all the log files of the form "G4_*.log"
# found in the same directory as the script.
# It looks at the name of the Worker Node (WN), which is
# found at the 4th line of the log file, and keep track of how
# many times the same WN has been used.
# Then t prints out the list of WN, ordered according to
# the number of times they have been used, and the corresponding
# number.
# Finally, we print also the "Site" list, i.e. the number of
# jobs that have been completed by each "site", where a "site"
# is defined here as what follows the first "." in the WN name:
#   e.g.  worker node name:  d0cs1101.fnal.gov
#         corresponding site name:    fnal.gov
# NB) The "site" is very close to "Computing Element" (CE),
#     but it is not identical. For instance:
#       osg-gw-2.t2.ucsd.edu:2119/jobmanager-condor-geant4
#       osg-gw-4.t2.ucsd.edu:2119/jobmanager-condor-geant4
#     are two different CEs, but they count as the same "site"
#     according to my definition, i.e. "t2.ucsd.edu".
#
#--------------------------------------------------------------------------

import os
import sys
import string
import math


#===============================================
#==================== MAIN ===================== 
#===============================================

# Check whether the number of arguments are correct or not.
if ( len( sys.argv ) >= 2 ) :
    print " Usage:  runningSites.py "
else :

    os.system( "ls -1 G4_*.log > listLogs.txt" )
    listLogs = open( "listLogs.txt" )
    dictWN = {}                                # empty dictionary
    dictSite = {}
    for fileName in listLogs :
        #print fileName.strip()                ###DEBUG
        thefile = open( fileName.strip() )
        counterLine = 0
        for line in thefile :
            #print line.strip()                ###DEBUG
            counterLine += 1

            if ( counterLine == 4 ) :          # the WN is at the 4th line

                name_wn = line.strip()
                #print name_wn                 ###DEBUG
                if dictWN.has_key( name_wn ) :
                    dictWN[name_wn] += 1
                else :
                    dictWN[name_wn] = 1

                name_site = name_wn[name_wn.find(".")+1:]
                #print 'name_site=', name_site ###DEBUG
                if dictSite.has_key(name_site) :
                    dictSite[name_site] += 1
                else :
                    dictSite[name_site] = 1

                break
        thefile.close()
    listLogs.close()
    os.system( "rm listLogs.txt" )

    print '=== WN list ===='

    #print ' dictWN=', dictWN                  ###DEBUG            

    listValues = dictWN.values()
    #print '  listValues=', listValues         ###DEBUG
    listValues.sort()
    #print '  listValues after sorting=', listValues   ###DEBUG
    listValues.reverse()
    #print '  listValues after sorting and reversing=', listValues ###DEBUG
    previous_num = 0
    for num in listValues :
        if num != previous_num :
            previous_num = num
            for k in dictWN.keys() :
                if num == dictWN[k] :
                    print num, k 

    print '  '
    print '=== SITE list ===='

    #print ' dictSite=', dictSite              ###DEBUG            

    listValues = dictSite.values()
    #print '  listValues=', listValues         ###DEBUG
    listValues.sort()
    #print '  listValues after sorting=', listValues   ###DEBUG
    listValues.reverse()
    #print '  listValues after sorting and reversing=', listValues ###DEBUG
    previous_num = 0
    for num in listValues :
        if num != previous_num :
            previous_num = num
            for k in dictSite.keys() :
                if num == dictSite[k] :
                    print num, k 


#-------------------------------------------------------------------
