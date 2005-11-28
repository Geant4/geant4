#!/usr/bin/python

#----------------------------------------------------------------
# Last update: 28-Nov-2005.  
#
# This Python script is used for post-processing analysis, i.e.
# to monitor the overall p-value distributions for all jobs.
# We are assuming here that the log files of these simulation
# runs are collected in a directory, eventually with a
# subdirectory structure. This script should be run the tree
# parent directory.
#
# This script does not have input arguments, but you have to set
# appropriately some parameters, that you can find below (search
# for the keyword  ***LOOKHERE*** ). These are the
# parameters you need to set:
#    1)  directory : where to look for the log files from which
#                    the information are extracted. It looks in
#                    the specified directory and, recursively,
#                    to all its subdirectories.
#    2)  tuplePhysicsLists : Geant4 Physics Lists to be considered.
#
# This script produces in output, in the same directory where this
# script is located, the following files:
#
#    o  listFoundPvaluesFiles.txt : list of found files
#                                   "outputPvalues.log-*"
#    o  filePvalues.txt-LHEP-AD   : lines representing the number
#                                   of entries for the 100 bins of
#                                   pvalues, between 0 and 1, for
#                                   the Anderson-Darling test. This
#                                   file is meant to be read by a Paw
#                                   kumac (as pvalues.kumac).
#    o  filePvalues.txt-LHEP-C2   : as above for the Chi2 test.
#    o  filePvalues.txt-LHEP-CVM  : as above for the Cramer-von Mises  test.
#    o  filePvalues.txt-LHEP-KS   : as above for Kolmogorov-Smirnov test.
#    o  listCreatedPvaluesFiles.txt : list of files created
#                                     by this script.
#
#----------------------------------------------------------------

import os
import sys
import string
import math
import random

#***LOOKHERE***
# Look for the files in this directory and recursively in all
# its subdirectories.
#directory = "/afs/cern.ch/sw/geant4/stat_testing/june05"
directory = "/users/ribon/dirAcceptanceSuite/dirObservables/dirTestPvaluesPY"

###tuplePhysicsLists = ("LHEP", "QGSP", "QGSC", "QGSP_BIC", "QGSP_BERT")
tuplePhysicsLists = ("LHEP",)

#***endLOOKHERE***

# Collect in files the following lists:
#  -  the list of the files created by this script;
#  -  the list of files that have been found;
createdFiles = open( "listCreatedPvaluesFiles.txt", "w" )
foundFiles = open( "listFoundPvaluesFiles.txt", "w" )

createdFiles.write( "listCreatedPvaluesFiles.txt" + "\n" )
createdFiles.write( "listFoundPvaluesFiles.txt" + "\n" )

#------------------------------------------------------------
# ---------------------   CLASSESS   ------------------------
#------------------------------------------------------------
class Distribution :
     # This class is the base class for any 1-dimensional distribution.
     # The vector of bin edges for the variable has to be specified
     # in the constructor.
     # It keeps count of the number of entries, the sum of the weights
     # for each bin, and the sum of the square of the weights for each
     # bin (which is needed for the calculation of the statistical error
     # in each bin), including underflows, and overflows.
     # It can be normalized, either uniformely for all bins, or
     # independently bin by bin. The normalization factor is
     # multiplied to the contents, and, by squaring it, to the
     # square sum of the weights, but not to the entries.

     def __init__( self, name_input = "Distribution", binEdgesVec_input = [] ) :
         # The vector of increasing bin edges must be specified in
         # the constructor.
         self.thename = name_input
         # Check that the values in binEdgesVec_input are increasing.
         previousX = -1.0e100   # A very large negative number
         isOK = 1
         for x in binEdgesVec_input :
             if ( previousX > x ) :
                 isOK = 0
             previousX = x
         self.binEdgesVec = []
         if ( not isOK ) :
             print " ***ERROR*** NOT INCREASING X BINS EDGES"
         else :
             self.binEdgesVec = binEdgesVec_input
         # The content vector has a length equal to the vector
         # of x bin edges plus 1, because the its first element [0]
         # will collect the underflow values, and the last one
         # will collect the overflow values; the rest is for the
         # content of the bins.
         # We keep also two similar vectors: one with the sum of
         # the square of the weight of each entry, which is useful
         # to calculate the statistical error; and another with
         # the number of entries (which coincides with the content
         # vector in the case of weights = 1 and not normalization
         # applied).
         self.contentVec = [ 0.0  for i in xrange( len( self.binEdgesVec ) + 1 ) ]
         self.weight2Vec = [ 0.0  for i in xrange( len( self.binEdgesVec ) + 1 ) ]
         self.entriesVec = [  0   for i in xrange( len( self.binEdgesVec ) + 1 ) ]
     
     def name( self ) :
         return self.thename

     def numberOfBins( self ) :
         return len( self.binEdgesVec ) - 1

     def fullBinEdgesVec( self ) :
         return self.binEdgesVec

     def onlyBinContentVec( self ) :
         return self.contentVec[1:-1]

     def onlyBinErrorVec( self ) :
         # The statistical error of each bin, in the more general case
         # of weighted distribution, is given by the square root of the
         # quadratic sum of the weight of each entry, which is kept in
         # the vector  weight2Vec . Notice that this also applied in
         # the simple case of  weight = 1 , in which it reduces simply
         # to the square root of the content of the bin (i.e. of the
         # number of entries).
         return [ math.sqrt( x ) for x in self.weight2Vec[1:-1] ]

     def onlyBinEntriesVec( self ) :
         return self.entriesVec[1:-1]

     def sum_onlyBinContentVec( self ) :
         sum = 0.0
         for y in self.contentVec[1:-1] :
             sum += y
         return sum

     def sum_onlyBinEntriesVec( self ) :
         sum = 0
         for y in self.entriesVec[1:-1] :
             sum += y
         return sum

     def fullContentVec( self ) :
         return self.contentVec

     def fullEntriesVec( self ) :
         return self.entriesVec

     def sum_fullContentVec( self ) :
         sum = 0.0
         for y in self.contentVec :
             sum += y
         return sum

     def sum_fullEntriesVec( self ) :
         sum = 0
         for y in self.entriesVec :
             sum += y
         return sum

     def underflowContent( self ) :
         return self.contentVec[0]

     def underflowEntries( self ) :
         return self.entriesVec[0]

     def overflowContent( self ) :
         return self.contentVec[ len( self.binEdgesVec ) ]

     def overflowEntries( self ) :
         return self.entriesVec[ len( self.binEdgesVec ) ]

     def reset( self ) :
         # Reset to zero content, sum weights, and entries.
         self.contentVec = [ 0.0  for i in  xrange( len( self.binEdgesVec ) + 1 ) ]
         self.weight2Vec = [ 0.0  for i in  xrange( len( self.binEdgesVec ) + 1 ) ]
         self.entriesVec = [ 0    for i in  xrange( len( self.binEdgesVec ) + 1 ) ]

     def normalize( self, factor = 1.0, binPosition = 696969 ) :
         # This method does the overall normalization, i.e. the
         # same scaling factor is applied to all bins, if binPosition
         # has its default value (chosen to be a very large and atypical
         # number) otherwise it normalizes only the specified bin (this
         # is sometimes needed when you have non-constant binnings).
         # The normalization affects the content, and the weight2,
         # but not the entries. Notice that the scaling of weight2
         # is such to preserve the relative statistical error (or
         # precision), i.e. the ratio between the statistical error
         # in a given bin and the content of the bin itself.
         if ( factor < 0.0 ) :
             print " ***ERROR*** NEGATIVE NORMALIZATION FACTOR ", factor
             return
         if (  binPosition != 696969 ) :  # Normalize only one bin
             if ( binPosition < 0  or
                  binPosition >= len( self.contentVec ) ) :
                 print " ***ERROR*** WRONG BIN TO NORMALIZE ", binPosition
             else :
                 print " Normalizing : binPosition=", binPosition
                 print "  before: content=", self.contentVec[ binPosition ], \
                       "  weight2=", self.weight2Vec[ binPosition ], \
                       "  entries=", self.entriesVec[ binPosition ]
                 self.contentVec[ binPosition ] *= factor
                 self.weight2Vec[ binPosition ] *= factor * factor
                 print "  after: content=", self.contentVec[ binPosition ], \
                       "  weight2=", self.weight2Vec[ binPosition ], \
                       "  entries=", self.entriesVec[ binPosition ]
         else :                           # Normalize all bins
             for i in xrange( len( self.contentVec ) ) :
                 print " Normalizing : i=", i
                 print "  before: content=", self.contentVec[ i ], \
                       "  weight2=", self.weight2Vec[ i ], \
                       "  entries=", self.entriesVec[ i ]
                 self.contentVec[ i ] *= factor
                 self.weight2Vec[ i ] *= factor * factor
                 print "  after: content=", self.contentVec[ i ], \
                       "  weight2=", self.weight2Vec[ i ], \
                       "  entries=", self.entriesVec[ i ]

     def fill( self, xval, weight = 1.0 ) :
         # This method fills the content of the bin corresponding
         # to the value  xval , with the specified  weight
         # (default is 1.0).
         # If the weight is negative or there are no bins available
         # then print an error message and exit.
         # If xval is less than the first x bin edge, then
         # it is an underflow (first bin of contentVec);
         # If xval is greater than the last x bin edge, then
         # it is an overflow (last bin of contentVec).
         if ( weight < 0.0  or  self.numberOfBins() <= 0 ) :
             if ( weight < 0.0 ) :
                 print " ***ERROR*** NEGATIVE WEIGHT ", weight
             else :
                 print " ***ERROR*** NO BINS "                 
             return
         index = 0
         for x in self.binEdgesVec :
             if ( xval < x ) :
                 self.contentVec[ index ] += weight
                 self.weight2Vec[ index ] += weight * weight
                 self.entriesVec[ index ] += 1
                 break
             index += 1
         else :
             self.contentVec[ index ] += weight
             self.weight2Vec[ index ] += weight * weight
             self.entriesVec[ index ] += 1

     def meanAndSigma( self ) :
         # This method calculates the  mean  and  sigma  of the
         # weighted distribution, and their errors.
         # The middle of each bin is considered for computing the mean,
         # whereas underflow and overflow entries are ignored. 
         sum = 0.0
         sum2 = 0.0
         sumW = 0.0
         for i in xrange( self.numberOfBins() ) :
             x = 0.5 * ( self.binEdgesVec[ i ] + self.binEdgesVec[ i + 1 ] )
             sum += x * self.contentVec[ i + 1 ]
             sum2 += x * x * self.contentVec[ i + 1 ]
             sumW += self.contentVec[ i + 1 ]
         ( mean, sigma_mean, sigma, sigma_sigma ) = ( 0.0, 0.0, 0.0, 0.0 )
         if ( sumW > 1.0 ) :
              mean = sum / sumW
              sigma = math.sqrt( ( sum2 - sumW*mean*mean ) / ( sumW - 1 ) )
              sigma_mean = sigma / math.sqrt( sumW )
              sigma_sigma = sigma / math.sqrt( 2 * ( sumW - 1 ) )
         else :
              print " WARNING in  meanAndSigma() : sumW <= 1 ", sumW
         return ( mean, sigma_mean, sigma, sigma_sigma )


class Pdistribution( Distribution ) :
    # This class inherits from Distribution, and represents the
    # p-value distribution.

    def __init__( self, name = "p-value" ) :
        # Constructor: 100 bins, of equal size (0.01), from 0.0 to 1.0 .
        self.binSize = 0.01
        self.p_bins = [ 0.0 + self.binSize*i  for i in xrange( 101 ) ]
        #print "p_bins=", self.p_bins
        self.numberOfNan = 0
        self.numberOfInf = 0
        Distribution.__init__( self, name, self.p_bins )

    def getNumberOfNan( self ) :
         return self.numberOfNan
    def increaseNumberOfNan( self ) :
         self.numberOfNan += 1
    
    def getNumberOfInf( self ) :
         return self.numberOfInf
    def increaseNumberOfInf( self ) :
         self.numberOfInf += 1
    

#-------------------------------------------------------------
# ---------------------   FUNCTIONS   ------------------------
#-------------------------------------------------------------

def printParameters() :
    # 
    # This function prints out the values of the various parameters.
    print '  --- Start function  printParameters  --- '
    print '  directory = ', directory
    print '  Physics Lists : '
    for iPL in tuplePhysicsLists :
        print '                  ' , iPL
    print '  --- End   function  printParameters  --- '
    return


def extractInfo( inputFile ) :
    # 
    # This function receives in input a file, and then it extracts
    # from it the distribution of p-values for 4 different
    # statistical tests:
    #   1) Chi2 : "C2"
    #   2) Kolmogorov-Smirnov : "KS"
    #   3) Cramer-von Mises : "CVM"
    #   4) Anderson-Darling : "AD"
    # The function returs these four distributions.
    # Notice that it takes into account also anomalous p-values,
    # like "nan" or "inf"; negative values, or values >= 1 are
    # also described as underflow and overflow entries.
       
    print '  --- Start function  extractInfo  --- '

    pValueDistC2  = Pdistribution( iPL + "-C2" )
    pValueDistKS  = Pdistribution( iPL + "-KS" )
    pValueDistCVM = Pdistribution( iPL + "-CVM" )
    pValueDistAD  = Pdistribution( iPL + "-AD" )

    # Just for debugging, fill the distributions with random numbers.
    #for n in xrange( 10000 ) :
    #     pValueDistC2.fill( random.random() ) 
    #     pValueDistKS.fill( random.random() ) 
    #     pValueDistCVM.fill( random.random() ) 
    #     pValueDistAD.fill( random.random() ) 

    for line in inputFile :
         #print line.strip()
         if ( line.find( "pvalue=" ) > -1 ) :
              testName = line.split()[0]
              pValueStr = line.split()[3].replace( 'pvalue=', '' )
              print " testName=", testName, " pValueStr=", pValueStr
              pValueDist = 0
              if ( testName.find( "Chi2" ) > - 1 ) :
                   pValueDist = pValueDistC2
              elif ( testName.find( "KS" ) > - 1 ) :
                   pValueDist = pValueDistKS
              elif ( testName.find( "CVM" ) > - 1 ) :
                   pValueDist = pValueDistCVM
              elif ( testName.find( "AD" ) > - 1 ) :
                   pValueDist = pValueDistAD
              else :
                   print " ***ERROR*** : not recognized testName=", testName
              if ( pValueDist ) :
                   print " pValueStr=", pValueStr
                   if ( pValueStr.find( "nan" ) > -1 ) :
                        pValueDist.increaseNumberOfNan()
                        print " Nan : ", pValueStr
                   elif ( pValueStr.find( "inf" ) > -1 ) :
                        pValueDist.increaseNumberOfInf()
                        print " Inf : ", pValueStr
                   else :
                        pValueDist.fill( float( pValueStr ) )
                        print " number : ", pValueStr

    print '  --- End   function  extractInfo  --- '
    
    return ( pValueDistC2, pValueDistKS, pValueDistCVM, pValueDistAD )


#-------------------------------------------------------------
# ---------------------   MAIN   -----------------------------
#-------------------------------------------------------------

print '  ========== START pvalues.py ========== '

printParameters();

for iPL in tuplePhysicsLists :
    print ' --- Physics List : ', iPL, ' --- ' 

    # Creates the text files to be read by Paw.
    filePvaluesC2 = open( "filePvalues.txt-" + iPL + "-C2", "w" ) 
    createdFiles.write( "filePvalues.txt-" + iPL + "-C2" + "\n" )
    filePvaluesKS = open( "filePvalues.txt-" + iPL + "-KS", "w" ) 
    createdFiles.write( "filePvalues.txt-" + iPL + "-KS" + "\n" )
    filePvaluesCVM = open( "filePvalues.txt-" + iPL + "-CVM", "w" ) 
    createdFiles.write( "filePvalues.txt-" + iPL + "-CVM" + "\n" )
    filePvaluesAD = open( "filePvalues.txt-" + iPL + "-AD", "w" ) 
    createdFiles.write( "filePvalues.txt-" + iPL + "-AD" + "\n" )

    # Look for the  outputPvalues.log-*  files in all the subdirectories
    # of the current directory.
    fileName = "outputPvalues.log-" + iPL.strip() + "-*" 
    command = "find " + directory + "/." + " -name " + fileName
    command = command + " -follow -print > thePath.txt"
    #print '  command = ', command
    os.system( command )

    if ( os.path.getsize( "thePath.txt" ) > 0 ) :
         # Read the full path of the input files from  thePath.txt
         # and then open the input file.
         thePathFile = open( "thePath.txt", "r" )
         fullNameInputFile = ""
         for line in thePathFile :
              print " file : ", line.strip()
              fullNameInputFile = line.strip()
              foundFiles.write( fullNameInputFile + "\n" )
              inputFile = open( fullNameInputFile, "r" )

              pValueDist_list = extractInfo( inputFile )
              for pValueDist in pValueDist_list :
                   # Print all the information on the distribution
                   print " --- Distribution : ", pValueDist.name(), " --- " 
                   print " onlyBinEntriesVec=", pValueDist.onlyBinEntriesVec()
                   print " numberOfNan=", pValueDist.getNumberOfNan()
                   print " numberOfInf=", pValueDist.getNumberOfInf()
                   print " underflowEntries=", pValueDist.underflowEntries()
                   print " overflowEntries=", pValueDist.overflowEntries()
                   print " numberOfBinEntries=", pValueDist.sum_onlyBinEntriesVec()
                   ( mean, rms_mean, rms, rms_rms ) = pValueDist.meanAndSigma()
                   print " mean=", mean, " +/- ", rms_mean, "  (expected:  0.5 )"
                   print " rms=", rms, " +/- ", rms_rms, "  (expected:  1/sqrt(12) = 0.288675 )"
                   if ( rms_mean > 0.0  and
                        math.fabs( mean - 0.5 ) / rms_mean > 3.0 ) :
                        print " MEAN more than  3 Sigma  away from the expected value "
                   if ( rms_rms > 0.0  and
                        math.fabs( rms - 1.0/math.sqrt( 12.0 ) ) / rms_rms > 3.0 ) :
                        print " RMS more than  3 Sigma  away from the expected value "                    
                   # Write the distribution on a text file to be read by Paw.
                   for val in pValueDist.onlyBinEntriesVec() :
                        if ( pValueDist.name().find( "C2" ) > - 1 ) :
                             filePvaluesC2.write( str( val ) + "\n" )
                        elif ( pValueDist.name().find( "KS" ) > - 1 ):
                             filePvaluesKS.write( str( val ) + "\n" )
                        elif ( pValueDist.name().find( "CVM" ) > - 1 ):
                             filePvaluesCVM.write( str( val ) + "\n" )
                        elif ( pValueDist.name().find( "AD" ) > - 1 ):
                             filePvaluesAD.write( str( val ) + "\n" )

    filePvaluesC2.close()
    filePvaluesKS.close()
    filePvaluesCVM.close()
    filePvaluesAD.close()
                 
# Close the files.
foundFiles.close()
createdFiles.close()
    
print '  ========== END pvalues.py ========== '


