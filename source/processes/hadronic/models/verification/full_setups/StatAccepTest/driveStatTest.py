#!/usr/bin/python

#----------------------------------------------------------------
# Last update: 13-Aug-2009.
#
# This Python script has three input parameters, which are
# strings that allows to infer the names of the two files
# containing ROOT trees. These trees are used by the Root
# macro mainRootScript.C , that does statistical tests
# between the various distributions.
#----------------------------------------------------------------

import os
import sys

print '    ========== BEGIN driveStatTest.py ========== '

if ( len( sys.argv ) > 3 ) :
    caseA = sys.argv[1]
    caseB = sys.argv[2]
    generalCase = sys.argv[3]
    if ( len( sys.argv ) > 4 ) :
        print '    TOO MANY ARGUMENTS: only the first three are considered! '

# Prepare the input files for the  pvalues  program.

file1 = "ntuple.root" + "-" + caseA + "-" + generalCase
print "    file1 = ", file1
if ( not os.path.exists( file1 ) ) :
    print '    ***ERROR*** in driveStatTest.py : file1=', file1, '  NOT found!'
    sys.exit( 51 )

file2 = "ntuple.root" + "-" + caseB + "-" + generalCase
print "    file2 = ", file2
if ( not os.path.exists( file2 ) ) :
    print '    ***ERROR*** in driveStatTest.py : file2=', file2, '  NOT found!'
    sys.exit( 52 )

os.system( "ln -sfn " + file1 + " ntuple_a.root" )
os.system( "ln -sfn " + file2 + " ntuple_b.root" )

fileLog = "outputPvalues.log" + "-" + generalCase   # Log file for pvalues.
print "    fileLog = ", fileLog

# Execute the  pvalue  executable that does the Statistical tests.
resultCode = os.system( "root -b -q mainRootScript.C > " + fileLog + " 2>&1 " )
if ( resultCode != 0 ) :
    print ' ***ERROR*** from: os.system( root -b -q mainRootScript.C ) ! code=', \
          resultCode
    sys.exit( 53 )

# Delete the symbolic links.
os.remove( "ntuple_a.root" )
os.remove( "ntuple_b.root" )

print '    ========== END driveStatTest.py ========== '
