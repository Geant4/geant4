#!/usr/bin/python

#----------------------------------------------------------------
# This Python script has one input parameter, a string,
# that represents the name of a file (in the same directory
# where this script is located) used as input in this script
# (the input file is output of the program  pvalue  that makes
#  several statistical tests between a certain number of
#  observables) to do the following.
# It reads the input file, and when it finds the string
# "***WARNING***"  it takes the name of the corresponding
# observable (that fails at least one statistical test) and
# finds the ID of the corresponding histogram that shows such
# an observable. This ID value is then passed to a  .kumac ,
# that is run in batch (not interactively!), and produces a
#  .ps  file with the plot of that observable.
# Finally, this file is renamed in such a way that the new name
# remembers which was the name of the input file and of the
# observable whose plot it represents.
#----------------------------------------------------------------

import os
import sys
 
print '     ========== BEGIN plot.py ========== '

inputFile = ""

if ( len( sys.argv ) > 1 ) :
    inputFile = sys.argv[1]
    if ( len( sys.argv ) > 2 ) :
        print '     TOO MANY ARGUMENTS: all but the first are ignored '

if ( not os.path.exists(inputFile) ) :
    print '     ***ERROR*** in plot.py : input file', inputFile, '  NOT found!'
    sys.exit(0)

pfile = open(inputFile, 'r')

# Extract from the input file the label to identify the case
# considered: it is the string after .log . For instance:
# outputPvalues.log-LHEP-FeSci-p-20GeV-5k  ->  -LHEP-FeSci-p-20GeV-5k
label = ""
if ( inputFile.find( ".log" ) ) :
    pos = inputFile.find( ".log" )
    label = inputFile[pos+4:]
###print '     label = ', label 

# --------------------------------------------------------
# Parse the specified file to find the failing comparisons 
# --------------------------------------------------------

for line in pfile :
###    print '     line = ', line
    if ( line.find ("***WARNING***") >= 0 ) :
###        print "     line.split() = ", line.split()
        numObservable = line.split()[1]
###        print '     numObservable =', numObservable

        # Schema for the id of the observables:
        #     Observable   1   ->   Histogram   11
        #     Observable   2   ->   Histogram   12
        #     Observable  L0   ->   Histogram  100
        #     Observable  L1   ->   Histogram  101
        #     ...
        #     Observable  L99  ->   Histogram  199
        #     Observable  R0   ->   Histogram  200
        #     Observable  R1   ->   Histogram  201
        #     ...
        #     Observable  R29  ->   Histogram  229
        
        id = ""
        if ( numObservable == "1" ) :
            id = "11"
        elif ( numObservable == "2" ) :
            id = "12"
        elif ( numObservable[0] == "L" ) :
            if ( len( numObservable ) > 2 ) :
                id = "1" + numObservable[1:3]
            else :
                id = "10" + numObservable[1]
        elif ( numObservable[0] == "R" ) :
            if ( len( numObservable ) > 2 ) :
                id = "2" + numObservable[1:3]
            else :
                id = "20" + numObservable[1]

        print '     id = ', id 
        if ( not id ) :
            print '     ***ERROR*** in plot.py : SOMETHING WRONG WITH THE ID : id = ', id
        else :

            #----------------------------------------------
            # Plot the observable that fails the comparison
            #----------------------------------------------

            # Create  main.kumac  that simply calls  plot.kumac
            # with the right argument: this is necessary to execute
            # paw in batch (no arguments are allowed).
            if ( os.path.exists("main.kumac") ) :
                os.system("rm main.kumac")
                
            mainKumac = open("main.kumac", 'w')
            command = "exec plot.kumac " + id + "\n"
            mainKumac.write(command)
            mainKumac.close();

            # Execute in batch the kumac (without arguments) main.kumac
            os.system("paw -w 0 -b main.kumac")

            # Rename the .ps file to remember the observable
            os.rename( "plot.ps" , "plot.ps-" + numObservable + label )

print '     ========== END plot.py ========== '
