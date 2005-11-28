#-------------------------------------------------------------------
# Last update: 28-Nov-2005
#
# This script should be run in the directory which has, as immediate
# subdirectories, the results of the Grid validation testing,
# which consist of tar-balls (one per subdirectory, each subdirectory
# being a job).
#
# The script does is following: it goes in each subdirectory, unpackes
# the tar-ball, and keep a statistics of the files that should be
# contained in such tar-ball. At the end, it prints out this
# information, and also produces the following files:
#    listDirFailed.txt : a list of the directories whose tar-ball
#                        do not have the expected files.
#    listPS.txt : a list of .ps files that should be examined.
#
# This script can be run also to get only the statistics, without
# unpacking the tar-balls: to do so, comment the line where
# ***LOOKHERE*** appears.
#
# This script has not input argument, and should be run as:
#         $  python unpack
#
#-------------------------------------------------------------------

import os
import sys
import string

listDirFailed = open( "listDirFailed.txt", "w" )
listPS = open( "listPS.txt", "w" )

os.system( "ls -1F | grep / > listDir.txt" )
listDir = open( "listDir.txt", "r" )

countDir = 0
countDirWithNtuples = 0
countDirWithPvalues = 0
countDirWithPS = 0
countPS = 0

for dir in listDir :
    #print dir
    countDir += 1 
    saveDir = os.getcwd()
    os.chdir( dir.strip() )
    currentDir = os.getcwd()
    #print " I am in directory: ", currentDir

    os.system( "tar xvfz *.tgz" )       #***LOOKHERE**

    os.system( "ls -1 > listFiles.txt" )
    listFiles = open( "listFiles.txt", "r" )
    foundNtuples = 0
    foundPvalues = 0
    foundPS = 0
    for iFile in listFiles :
        #print " iFile=", iFile
        if ( iFile.find( "ntuple.hbook" ) > -1 ) :
            foundNtuples = 1
        if ( iFile.find( "outputPvalues.log" ) > -1 ) :
            foundPvalues = 1
        if ( iFile.find( "plot.ps" ) > -1 ) :
            foundPS = 1
            countPS += 1
            listPS.write( currentDir + "/" + iFile )
    if ( foundNtuples ) :
        countDirWithNtuples += 1
    if ( foundPvalues ) :
        countDirWithPvalues += 1
    if ( foundPS ) :
        countDirWithPS += 1
    if ( not foundNtuples  or  not foundPvalues ) :
        listDirFailed.write( currentDir + "\n" )

    listFiles.close()
    os.system( "rm listFiles.txt" )
    os.chdir( saveDir)

print " Summary: "
print " Number of directories = ", countDir
print " Number of directories with Ntuples = ", countDirWithNtuples
print " Number of directories with p-values = ", countDirWithPvalues
print " Number of directories with .PS = ", countDirWithPS
print " Number of .PS = ", countPS

listDir.close()
listDirFailed.close()
listPS.close()


#-------------------------------------------------------------------

