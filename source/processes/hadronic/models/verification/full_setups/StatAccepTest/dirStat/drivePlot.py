#!/usr/bin/python

#----------------------------------------------------------------
# This Python script has three input parameters, which are
# strings that allows to infer the names of the two files
# containing HBOOK ntuples. These ntuples are used by the
# program  pvalues , that does statistical tests between
# the various distributions.
# This script  drivePlot.py  drives the execution of the
# program  pvalues , and then drives the execution of another
# Python script,  plot.py , passing to it the name of the log
# file of  pvalues .
#----------------------------------------------------------------

import os
import sys

if ( len( sys.argv ) > 3 ) :
    caseA = sys.argv[1]
    caseB = sys.argv[2]
    generalCase = sys.argv[3]
    if ( len( sys.argv ) > 4 ) :
        print ' TOO MANY ARGUMENTS: only the first three are considered! '

# Prepare the input files for the  pvalues  program.

file1 = "../ntuple.hbook" + "-" + caseA + "-" + generalCase
print "file1 = ", file1
if ( not os.path.exists( file1 ) ) :
    print '***ERROR*** : file1=', file1, '  NOT found!'
    sys.exit(0)

file2 = "../ntuple.hbook" + "-" + caseB + "-" + generalCase
print "file2 = ", file2
if ( not os.path.exists( file2 ) ) :
    print '***ERROR*** : file2=', file2, '  NOT found!'
    sys.exit(0)

os.system( "ln -s " + file1 + " ntuple_a.hbook" )
os.system( "ln -s " + file2 + " ntuple_b.hbook" )

fileLog = "outputPvalues.log" + "-" + generalCase   # Log file for pvalues.
print "fileLog = ", fileLog

# Execute the  pvalue  executable that does the Statistical tests.
os.system( "pvalue > " + fileLog + " 2>&1 " )

# Execute the Python scrip  plot.py  which uses the log file
# of the previous executable.
os.system( "python2.2 plot.py " + fileLog )

# Delete the symbolic links, and renamed the output files
os.remove( "ntuple_a.hbook" )
os.remove( "ntuple_b.hbook" )

os.rename( "histo_a.hbook", "histo_a.hbook" + "-" + generalCase )
os.rename( "histo_b.hbook", "histo_b.hbook" + "-" + generalCase )

os.rename( "cloudsA.xml", "cloudsA.xml" + "-" + generalCase )
os.rename( "cloudsB.xml", "cloudsB.xml" + "-" + generalCase )
