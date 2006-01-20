//===================================================================
// Last update: 20-Jan-2006                 
//
//                   Program  power.cpp
//                   ------------------
// 
// This program computes the power of one or more statistical tests
// with respect to a pair of parent distributions:
//                   ( parent1 , parent2 ) 
// and for few confidence levels.
// 
// The statistical tests which you can select are: 
//      1)  KolmogorovSmirnov
//      2)  CramerVonMisesUnbinned
//      3)  AndersonDarlingUnbinned
//      4)  Chi2
//      5)  CramerVonMisesBinned
//      6)  AndersonDarlingBinned
//      7)  Combined one: the best between 1), 2), 3), 4), as
//                        it is currently done in the Acceptance Suite. 
//
// The parent distributions that can be used are:
//      0)  VisibleEnergy :                  for the parent1;          
//      1)  VisibleEnergy_shifted( alpha ) : for the parent2,
//                                           where alpha is a parameter
//                                             between 0 and 1.
// See below for more details on these distributions.
// Notice that if you want to use this program to evaluate the statistical
// power for a different observable and/or a different use case, then it 
// should be enough to replace the two parent distributions "VisibleEnergy"
// and "VisibleEnergy_shifted( alpha )" with your preferred "MyObservable"
// and "MyObservable_shifted( alpha )".
//
// The confidence levels that are used by default are:
//      1)  90%
//      2)  95%
//      3)  99%
//      4)  99.9%
//
// The method is based on Monte Carlo pseudoexperiments, each 
// consisting in drawing two samples, S1 and S2, of equal size 
// from the parent distributions  parent1  and  parent2 , 
// respectively.
// Given then the confidence_level, the power of a given 
// statistics test is defined as the fraction of pseudoexperiments
// in which the p-value for the pair of samples (S1, S2) is
// less than  1 - confidence_level , i.e. :
// 
//              # pseudoexperiments with p-value < 1 - confidence_level 
//  power =  -----------------------------------------------------------
//                  # pseudoexperiments
 // 
// The power is computed here using only the analytical p-value 
// of the statistics test, usually based on some asympthotic 
// approximations. 
//
// You can select the various options by looking to string
// " ***LOOKHERE*** " in the program below.
//  
// Notice that, in practice, the only major parameter which
// needs to be changed explicitly by hand, and therefore 
// requires to rebuild the program and to re-execute it, is
// the sample size (which is the same for the two samples),
// and the shifting parameter (alpha) for the visible energy.
// 
// The output of the program consists of the printing on the 
// screen.
//
// NB) This program is a modified (simplified) version of the
//     power.cpp  program distributed with the  StatisticalToolkit.
// 
//===================================================================

#include <memory>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <cmath>
#include <string>

#include "AIDA/AIDA.h"

#include "KolmogorovSmirnovComparisonAlgorithm.h"
#include "Chi2ComparisonAlgorithm.h"
#include "AndersonDarlingUnbinnedComparisonAlgorithm.h"
#include "AndersonDarlingBinnedComparisonAlgorithm.h"
#include "CramerVonMisesUnbinnedComparisonAlgorithm.h"
#include "CramerVonMisesBinnedComparisonAlgorithm.h"
#include "ComparisonResult.h"
#include "StatisticsTesting/StatisticsComparator.h"

#include "CLHEP/Random/RanluxEngine.h"
#include "CLHEP/Random/RandGaussT.h"
#include "CLHEP/Random/RandBreitWigner.h"

using namespace StatisticsTesting; 


//==================================================================================
// CLASSES
//==================================================================================

// -----------------  Class  DistributionGenerator  ----------------------

// Class responsible for generating the two distributions.
// Notice that if you want to use this program to evaluate the statistical
// power for a different observable and/or a different use case, then it 
// should be enough to replace the two parent distributions "VisibleEnergy"
// and "VisibleEnergy_shifted( alpha )" with your preferred "MyObservable"
// and "MyObservable_shifted( alpha )".

class DistributionGenerator {

public:

  enum { NumberOfDistributions = 2 };

  enum Type { VisibleEnergy, VisibleEnergy_shifted };

  DistributionGenerator( const Type typeIn = VisibleEnergy );

  ~DistributionGenerator();

  double draw() const;
  // It draws a random number according to the selected distribution.

  static double getAlpha() { return alpha; }
  // Return the parameter alpha.

private:

  DistributionGenerator( const DistributionGenerator & x );
  // No copy constructor.  

  Type theType;

  static const double alpha, cA, sA, xjA, cB, sB, xjB, c, s, xj; 
  static const double pAVec[], pBVec[], pVec[];
  // Parameters needed for the visible energy distribution.
};

//***LOOKHERE***
const double DistributionGenerator::alpha = 0.5;  //***MAIN-PARAMETER***


// Constants for the visible energy distributions.
const double DistributionGenerator::cA = 8.9573;
const double DistributionGenerator::sA = -0.48792E-01;
const double DistributionGenerator::xjA = 70.0;

const double DistributionGenerator::cB = 9.1410;
const double DistributionGenerator::sB = -0.46462E-01;
const double DistributionGenerator::xjB = 75.0;

const double DistributionGenerator::c = 
  (1.0 - DistributionGenerator::alpha) * DistributionGenerator::cA  + 
  DistributionGenerator::alpha * DistributionGenerator::cB; 

const double DistributionGenerator::s = 
  (1.0 - DistributionGenerator::alpha) * DistributionGenerator::sA  + 
  DistributionGenerator::alpha * DistributionGenerator::sB; 

const double DistributionGenerator::xj = 
  (1.0 - DistributionGenerator::alpha) * DistributionGenerator::xjA  + 
  DistributionGenerator::alpha * DistributionGenerator::xjB; 

const double DistributionGenerator::pAVec[ 6 ] = {
  15.291 , 
  -15.217 , 
  4.3102 ,
  -0.14385 ,  
  0.17296E-02 ,
  -0.71385E-05 
};

const double DistributionGenerator::pBVec[ 6 ] = {
  10.942 ,
  -9.4891 ,
  2.2055 ,
  -0.51487E-01 , 
  0.34915E-03 ,
  -0.31766E-06 
};

const double DistributionGenerator::pVec[ 6 ] = {
  (1.0 - DistributionGenerator::alpha) * DistributionGenerator::pAVec[0] +
  DistributionGenerator::alpha * DistributionGenerator::pBVec[0] ,
  (1.0 - DistributionGenerator::alpha) * DistributionGenerator::pAVec[1] +
  DistributionGenerator::alpha * DistributionGenerator::pBVec[1] ,
  (1.0 - DistributionGenerator::alpha) * DistributionGenerator::pAVec[2] +
  DistributionGenerator::alpha * DistributionGenerator::pBVec[2] ,
  (1.0 - DistributionGenerator::alpha) * DistributionGenerator::pAVec[3] +
  DistributionGenerator::alpha * DistributionGenerator::pBVec[3] ,
  (1.0 - DistributionGenerator::alpha) * DistributionGenerator::pAVec[4] +
  DistributionGenerator::alpha * DistributionGenerator::pBVec[4] ,
  (1.0 - DistributionGenerator::alpha) * DistributionGenerator::pAVec[5] +
  DistributionGenerator::alpha * DistributionGenerator::pBVec[5] 
};


DistributionGenerator::DistributionGenerator( const Type typeIn ) : theType( typeIn ) {}


DistributionGenerator::~DistributionGenerator() {}


double DistributionGenerator::draw() const {

  double value = 0.0;

  switch ( theType ) {
  case VisibleEnergy :
    {
      // The distribution  VisibleEnergy  corresponds to a parametrization 
      // of the visible energy distribution in a Pb-Sci calorimeter, for a 
      // 9GeV pi+, obtained with Geant4 8.0, Physics List QGSP_GN, with 
      // default 0.7 mm production cut, running 20,000 simulated events. 
      // This distribution has been parametrized as the following splinline:
      //               | 0                           for   0 <= x < 2.0 MeV ;
      // function(x) = | polynomial of 5th degrees , for   2 <= x < xjA MeV ;
      //               | exp(c+s*x)                  fir   xjA <= x < 200.0 MeV.
      // The x variable is in MeV.
      // To generate values according to this function, whose maximum is 
      // (a bit) less than 900, we do the following: a candidate x value 
      // is draw from a flat random number between 0 and 200 (MeV); then 
      // we compute the corresponding value of the function (at this candidate
      // x value); then we draw another flat random number between 0 and 1 and 
      // compare it with the function value divided by 900 : if the 
      // random number is below such value we accept the candidate x value,
      // otherwise we repeat the procedure with another candidate x value.
      // If you want to plot the values that are generated, you have to
      // uncomment the line with the printing of x; then call output.log
      // the output produced by running this program; delete the lines
      // which are not values of x; and finally use the kumac  check.kumac .

      double x, r, fx;
      do {
      	x = 200.0 * HepRandom::getTheEngine()->flat();
      	if ( x < 2.0 ) {
      	  fx = 0.0;
      	} else if ( x < xjA ) {
      	  fx = pAVec[0] + pAVec[1]*x + pAVec[2]*x*x + pAVec[3]*x*x*x
      	     + pAVec[4]*x*x*x*x + pAVec[5]*x*x*x*x*x;
      	} else {
      	  fx = exp( cA + sA*x );
      	}
      	//std::cout << "      trying: x=" << x << "  f(x)=" << fx << std::endl;
      } while ( HepRandom::getTheEngine()->flat() > fx/900.0 );
      value = x;
      //std::cout << x << std::endl; //***DEBUG***
      break;
    }

  case VisibleEnergy_shifted :
    {
      // The family of distributions  VisibleEnergy_shifted( alpha ) ,  
      // where  alpha  varies between 0 and 1, corresponds to splinlines
      // as described above, and such that:
      //   VisibleEnergy_shifted( 0 ) = VisibleEnergy distribution
      //                                (see definition above);
      //   VisibleEneryg_shifted( 1 ) = parametrization of the visible energy 
      //                                distribution obtained for the same 
      //   configuration (9 GeV pi+ in PbSci, QGSP_GN, 0.7mm cut, 20k events), 
      //   but with the hadronic elastic G4LElasticB.cc .
      // For a given alpha, we use the same method as described above to
      // generate values according to this function, given the fact that
      // for all alpha values the maximum of the function is always (a bit) 
      // less than 900.
      // If you want to plot the values that are generated, you have to
      // uncomment the line with the printing of x; then call output.log
      // the output produced by running this program; delete the lines
      // which are not values of x; and finally use the kumac  check.kumac .

      double x, r, fx;
      do {
	x = 200.0 * HepRandom::getTheEngine()->flat();
	if ( x < 2.0 ) {
	  fx = 0.0;
	} else if ( x < xj ) {
	  fx = pVec[0] + pVec[1]*x + pVec[2]*x*x + pVec[3]*x*x*x
	     + pVec[4]*x*x*x*x + pVec[5]*x*x*x*x*x;
	} else {
	  fx = exp( c + s*x );
	}
	//std::cout << "\t trying: x=" << x << "  f(x)=" << fx << std::endl;
      } while ( HepRandom::getTheEngine()->flat() > fx/900.0 );
      value = x;
      //std::cout << "\t" << x << std::endl; //***DEBUG***
      break;
    }
  default : 
    {
      std::cout << "\t ***WRONG*** I should NOT be inside the  default: case" 
		<< std::endl; 
      break;
    }
  }

  return value;
}


// -----------------  Class  HistoBinningSelector  ----------------------------

// This class represents the binning of a histogram.

class HistoBinningSelector {

public:

  HistoBinningSelector();

  ~HistoBinningSelector();

  void easiest( const int n, const double inf, const double sup,
		std::vector< double > & binsVector );
  // The easiest binning, namely  n  bins of equal size between 
  //  inf  and  sup , is returned.
 
private:

  HistoBinningSelector( const HistoBinningSelector & x );
  // No copy constructor.  

  HistoBinningSelector & operator=( const HistoBinningSelector & x );  
  // No assignment operator.

};


HistoBinningSelector::HistoBinningSelector() {
  //std::cout << " Constructor of  HistoBinningSelector " << std::endl;
}


HistoBinningSelector::~HistoBinningSelector() {}


void HistoBinningSelector::easiest( const int n, const double inf, const double sup, 
				    std::vector< double > & binsVector ) {
  binsVector.clear();
  if ( n > 0  &&  inf < sup ) {
    double step = ( sup - inf ) / static_cast< double >( n ); 
    double x = inf;
    for ( int iBin = 0; iBin <= n; iBin++ ) {
      binsVector.push_back( x );
      x += step;
    }
  }
}


// ----------  Class GeneratedPseudoExpStruct ----------

// This is a struct which contains the pseudoexperiment number,
// the id (a number) identifing the parent distributions of the 
// first sample and of the second one (these two parent distributions 
// can be different, but always with the same sample size), the 
// generated values of the two samples, and the binning (which will 
// be used only for binned statistical tests).

struct GeneratedPseudoExpStruct {
  int pseudoExpNum;
  int idParents;
  std::vector< double > firstSample;
  std::vector< double > secondSample;
  std::vector< double > binning;
};


// ----------  Class StatisticsComparatorStruct ----------

// This is a struct which contains the information about the
// active statistics comparator: its id (a number), the full
// name (useful for printing out the results at the end),
// and whether it is binned or not.

struct StatisticsComparatorStruct {
  int id;
  std::string fullName;
  bool isBinned;
};


// ----------  Class TestResultPseudoExpStruct ----------

// This is a struct which contains the result of a statistical
// test on two samples (of equal size): this result consists
// of the pseudoexperiment number, an identifier of the statistical 
// test that has been considered, and the distance and the p-value 
// returned by that statistical test.

struct TestResultPseudoExpStruct {
  int pseudoExpNum;
  int idParents;
  int idStatTest;
  double distance;
  double pValue;
  double ndf;
};


// ----------  Class  PowerCalculator  ----------

// This class calculates the power of one or more statistical tests 
// by generating pseudoexperiments, in each of which the two samples 
// are drawn from two different parent distributions, and then the power 
// is the fraction of time the p-value is below a certain acceptance level, 
// i.e. the fraction of time a test recognizes them as coming from 
// different distributions.

class PowerCalculator {
  
public:

  PowerCalculator();
  ~PowerCalculator();

  inline bool isBinnedTest() const;
  inline void setIsBinnedTest( const bool isBinnedTestIn );
  // Get/Set the bool that tells you whether the statistical test 
  // we are examing is or not a binned test. 

  inline int currentPseudoExp() const;
  // It returns the "ordering number", from  0  to  PowerCalculator::numberPseudoExps-1, 
  // of the current generation of pair of distributions to be compared.
  
  void generateOnePseudoExp( std::vector< GeneratedPseudoExpStruct > & generation );
  // It fills the initial empty vector  generation  with all the considered
  // combinations of the first and second parent distributions. For each
  // of them, there are the two samples of equal size, and the binning
  // (which will be used only by binned statistical tests).

  void storeTestResultPseudoExp( const std::vector< TestResultPseudoExpStruct > & 
				 resultIn );
  // After that the previous method hands over the generated data to
  // "external" classes which do the actual statistical test, it needs
  // back the results of the various statistics tests, in order to calculate, 
  // by the next method, the power of the statistical test. 
  // Because we are considering the possibility of performing different 
  // statistical tests on the same two samples, a vector is the argument.

  void printResults( const std::vector< StatisticsComparatorStruct > & comparatorVec );
  // It prints out the power of the tests.

  static const int numberPseudoExps; // Number of Monte Carlo pseudoexperiments.

private:

  PowerCalculator( const PowerCalculator & x );
  // No copy constructor.  

  PowerCalculator & operator=( const PowerCalculator & x );  
  // No assignment operator.

  int theCurrentPseudoExp;
  bool theIsFailed;
  bool theIsBinnedTest;

  std::vector< TestResultPseudoExpStruct > collecResults;

  std::vector< bool > isOnParentDistribution;
 
  static const bool isOnVisibleEnergyParent;       
  static const bool isOnVisibleEnergyShiftedParent;

  static const double confidenceLevel1;       // Four confidence levels are
  static const double confidenceLevel2;       // possible, but you can eventually
  static const double confidenceLevel3;       // add other ones.
  static const double confidenceLevel4;

  static const int sampleSize;  // Number of elements (random drawings) in each sample.
  static const int sampleSizeParent;  // Number of elements (random drawings) for the
                                      // reference parent higher statistics sample.

  // Below the static constants relevant only for binned tests. 
  static const int numberBinsHisto;      // Number (max) of bins.
  static const double infHisto;          // Lowest value in the case of constant binning.
  static const double supHisto;          // Highest value "  "   "    "    "       "

  HistoBinningSelector theHistoBinningSelector;

};


//***LOOKHERE***
const bool PowerCalculator::isOnVisibleEnergyParent  = true;
const bool PowerCalculator::isOnVisibleEnergyShiftedParent  = true;

const double PowerCalculator::confidenceLevel1 = 0.90; 
const double PowerCalculator::confidenceLevel2 = 0.95; 
const double PowerCalculator::confidenceLevel3 = 0.99; 
const double PowerCalculator::confidenceLevel4 = 0.999; 

const int PowerCalculator::sampleSize       = 3000;   //***MAIN-PARAMETER***
const int PowerCalculator::numberPseudoExps = 1000;

const int PowerCalculator::numberBinsHisto = 100;
const double PowerCalculator::infHisto     = 0.0;
const double PowerCalculator::supHisto     = 200.0;


PowerCalculator::PowerCalculator() :  
  theCurrentPseudoExp( 0 ), theIsBinnedTest( false ) {

  isOnParentDistribution.push_back( isOnVisibleEnergyParent );
  isOnParentDistribution.push_back( isOnVisibleEnergyShiftedParent );

  std::cout << " Constructor of  PowerCalculator " << std::endl;
  for ( int i = 0; i < isOnParentDistribution.size(); i++ ) {
    std::cout << "\t isOnParentDistribution[" << i << "] : " 
	      << isOnParentDistribution[ i ] << std::endl;          
  }
  std::cout << "\t PowerCalculator::confidenceLevel1 = " 
	    << PowerCalculator::confidenceLevel1 << std::endl
            << "\t PowerCalculator::confidenceLevel2 = " 
	    << PowerCalculator::confidenceLevel2 << std::endl
            << "\t PowerCalculator::confidenceLevel3 = " 
	    << PowerCalculator::confidenceLevel3 << std::endl
            << "\t PowerCalculator::confidenceLevel4 = " 
	    << PowerCalculator::confidenceLevel4 << std::endl
            << "\t PowerCalculator::sampleSize = " 
	    << PowerCalculator::sampleSize << "\t <--- " << std::endl
	    << "\t PowerCalculator::numberPseudoExps = " 
	    << PowerCalculator::numberPseudoExps << std::endl
            << "\t Only for binned tests: " << std::endl
            << "\t \t PowerCalculator::numberBinsHisto = "
            << PowerCalculator::numberBinsHisto << std::endl
            << "\t \t PowerCalculator::infHisto = "
            << PowerCalculator::infHisto << std::endl 
            << "\t \t PowerCalculator::supHisto = "
            << PowerCalculator::supHisto << std::endl; 

  std::cout << " DistributionGenerator infos : "
            << "\t NumberOfDistributions = " 
	    << DistributionGenerator::NumberOfDistributions << std::endl
            << "\t alpha   (parametrization parameter) = " 
	    << DistributionGenerator::getAlpha() << "\t <--- " << std::endl;
}


PowerCalculator::~PowerCalculator() {} 


inline int PowerCalculator::currentPseudoExp() const {
  return theCurrentPseudoExp;
}


inline bool PowerCalculator::isBinnedTest() const {
  return theIsBinnedTest;
}

inline void PowerCalculator::setIsBinnedTest( const bool isBinnedTestIn ) {
  theIsBinnedTest = isBinnedTestIn;
}


void PowerCalculator::
storeTestResultPseudoExp( const std::vector< TestResultPseudoExpStruct > & resultIn ) {
  for ( std::vector< TestResultPseudoExpStruct >::const_iterator it = resultIn.begin();
	it != resultIn.end(); ++it ) {
    collecResults.push_back( *it );
  }
}


void PowerCalculator::
generateOnePseudoExp( std::vector< GeneratedPseudoExpStruct > & generation ) {

  if ( currentPseudoExp() == 0 ) {
    std::cout << " Generating pseudoexperiments..." << std::endl;
  }

  if ( currentPseudoExp() >= PowerCalculator::numberPseudoExps ) return;

  int n = currentPseudoExp(); 

  std::vector< std::vector< double > > firstGenerationVec;

  for ( int parentDist = 0; 
	parentDist < DistributionGenerator::NumberOfDistributions; parentDist++ ) {
    std::vector< double > firstSample;
    std::vector< double > secondSample;
    if ( PowerCalculator::isOnParentDistribution[ parentDist ] ) {
      DistributionGenerator generator;
      switch ( parentDist ) {
      case 0 : 
	generator = DistributionGenerator( DistributionGenerator::VisibleEnergy ); 
	break;
      case 1 : 
	generator = 
	  DistributionGenerator( DistributionGenerator::VisibleEnergy_shifted ); 
	break;
      default : 
	std::cout << "\t ***WRONG*** I should NOT be inside the  default: case" 
		<< std::endl; break;
      }
    
      // Draw the two samples with the same sample size. 
      for ( int i=0; i < PowerCalculator::sampleSize; i++ ) {
	firstSample.push_back( generator.draw() ); 
      }
    }
    firstGenerationVec.push_back( firstSample );
  }

  // Loop over the first parent distributions.
  for ( int parentDist1 = 0; 
	parentDist1 < DistributionGenerator::NumberOfDistributions; parentDist1++ ) {
    if ( ! PowerCalculator::isOnParentDistribution[ parentDist1 ] ) continue;

    // Loop over the second parent distributions.
    for ( int parentDist2 = parentDist1+1; 
	  parentDist2 < DistributionGenerator::NumberOfDistributions; parentDist2++ ) {
      if ( ! PowerCalculator::isOnParentDistribution[ parentDist2 ] ) continue;

      GeneratedPseudoExpStruct theGenStruct;
      theGenStruct.pseudoExpNum = n;
      theGenStruct.idParents    = parentDist1 + 100*parentDist2;
      theGenStruct.firstSample  = firstGenerationVec[ parentDist1 ]; 
      theGenStruct.secondSample = firstGenerationVec[ parentDist2 ]; //ALB-18Jan06 

      // Now determine the binning, if necessary.
      if ( isBinnedTest() ) {
	theHistoBinningSelector.easiest( PowerCalculator::numberBinsHisto,
					 PowerCalculator::infHisto,
					 PowerCalculator::supHisto,
					 theGenStruct.binning );
      }

      generation.push_back( theGenStruct );  

    } // End loop on the second parent
  } // End loop on the first parent

  // Print out a graphical aid to see how long the generation takes.
  if ( n % PowerCalculator::numberPseudoExps  == 0 ) {
    if ( PowerCalculator::numberPseudoExps >= 1000 ) {
      std::cout << "0%                                            100% " 
		<< std::endl;
    }
  }
  if ( PowerCalculator::numberPseudoExps >= 1000  &&
       ( n + 1 ) % ( PowerCalculator::numberPseudoExps / 50 )  == 0 ) {
    std::cout << "=" << std::flush;
    if ( ( n + 1 ) % PowerCalculator::numberPseudoExps  ==  0 ) {
      std::cout << std::endl;
    }
  }

  theCurrentPseudoExp++;
}


void PowerCalculator::
printResults( const std::vector< StatisticsComparatorStruct > & comparatorVec ) {

  // Loop over all possible pairs of parent distributions.
  for ( int parentDist1 = 0; 
	parentDist1 < DistributionGenerator::NumberOfDistributions; parentDist1++ ) {
    if ( ! PowerCalculator::isOnParentDistribution[ parentDist1 ] ) continue;
    for ( int parentDist2 = parentDist1+1; 
	  parentDist2 < DistributionGenerator::NumberOfDistributions; parentDist2++ ) {
      if ( ! PowerCalculator::isOnParentDistribution[ parentDist2 ] ) continue;

      std::string name1 = "Unknown";
      switch ( parentDist1 ) {
      case 0 : name1 = "VisibleEnergy"; break;
      case 1 : name1 = "VisibleEnergy_shifted"; break;
      }
      std::string name2 = "Unknown";
      switch ( parentDist2 ) {
      case 0 : name2 = "VisibleEnergy"; break;
      case 1 : name2 = "VisibleEnergy_shifted"; break;
      }
      std::cout << " ============================================= " << std::endl
		<< " Parent1 = "    << name1 << "\t ; Parent2 = " << name2 << std::endl;

      // Loop over the statistics tests.
      for ( std::vector< StatisticsComparatorStruct >::const_iterator itStatTest =
	      comparatorVec.begin(); itStatTest != comparatorVec.end(); ++itStatTest ) {

	std::cout << " --------------------------------------------- " << std::endl
		  << itStatTest->fullName << std::endl;

        // Loop over the collection of results, and fill the vector
        // theDistanceAndPValueCollection  only for those cases with
        // the current parent distributions for the two samples, and
        // with the current statistic test.
        std::vector< std::pair< double, double > > theDistanceAndPValueCollection;
	for ( std::vector< TestResultPseudoExpStruct >::const_iterator itResults =
		collecResults.begin(); itResults != collecResults.end(); ++itResults ) { 

	  if ( ( itResults->idParents == (parentDist1 + 100*parentDist2) ) &&
               ( itResults->idStatTest == itStatTest->id ) ) { 
	    double distance = itResults->distance;
	    double pValue   = itResults->pValue;
	    double ndf      = itResults->ndf;
	    theDistanceAndPValueCollection.push_back( std::pair< double, double >
						      ( distance,  pValue ) );
	  }
	}

	int numBelowConfidenceLevel1P = 0;
	int numBelowConfidenceLevel2P = 0;
	int numBelowConfidenceLevel3P = 0;
	int numBelowConfidenceLevel4P = 0;

	// Loop over the collection of pairs of distances and p-values.
	std::vector< std::pair< double, double > >::const_iterator iterCollection = 
	  theDistanceAndPValueCollection.begin();
	for ( int iPseudoExp = 0; 
	      iPseudoExp < PowerCalculator::numberPseudoExps; iPseudoExp++ ) {
	  
	  double p = iterCollection->second;
	  if ( p < -1.0E-9  ||  p - 1.0 > 1.0E-9 ) {
	    theIsFailed = true;
	    std::cout << "\t ***WRONG*** P-VALUES NOT INSIDE (0,1) " << p << std::endl; 
	  }  
	  
	  // Count the number of pseudoexperiments with p-values below 
          // 1 - confidenceLevel, for each of the four levels considered.
	  if ( p < (1.0 - PowerCalculator::confidenceLevel1) ) {
	    numBelowConfidenceLevel1P++;
	  }
	  if ( p < (1.0 - PowerCalculator::confidenceLevel2) ) {
	    numBelowConfidenceLevel2P++;
	  }
	  if ( p < (1.0 - PowerCalculator::confidenceLevel3) ) {
	    numBelowConfidenceLevel3P++;
	  }
	  if ( p < (1.0 - PowerCalculator::confidenceLevel4) ) {
	    numBelowConfidenceLevel4P++;
	  }

	  ++iterCollection;
	} // End loop over pseudoexperiments 
	
        double Nd = 1.0 * PowerCalculator::numberPseudoExps;
        // Calculate finally the power of the test for the four
        // confidence levels considered.
        for ( int level = 1; level <= 4; level++ ) { 
	  float cd, nd;
	  if ( level == 1 ) {
	    cd = PowerCalculator::confidenceLevel1;
	    nd = 1.0 * numBelowConfidenceLevel1P;
	  } else if ( level == 2 ) {
	    cd = PowerCalculator::confidenceLevel2;
	    nd = 1.0 * numBelowConfidenceLevel2P;
	  } else if ( level == 3 ) {
	    cd = PowerCalculator::confidenceLevel3;
	    nd = 1.0 * numBelowConfidenceLevel3P;
	  } else {
	    cd = PowerCalculator::confidenceLevel4;
	    nd = 1.0 * numBelowConfidenceLevel4P;
	  }
	  std::cout << "\t Confidence level=" << 100*cd << " % " << std::endl
                    << "\t \t \t powerP   = "  << 100*nd/Nd 
		    << "  +/-  " << 100 * sqrt( nd*(Nd-nd) / (Nd*Nd*Nd) ) 
                    << " % " << std::endl;
	}	
      } // End loop over the statistical tests that have been switched on
    } // End loop over the second parent distribution
  } // End loop over the first parent distribution

}


//==================================================================================
// GLOBAL FUNCTIONS
//==================================================================================

//***LOOKHERE*** Choose one or more statistical test to perform.
const bool isOnKolmogorovSmirnov       = true;
const bool isOnCramerVonMisesUnbinned  = true;
const bool isOnAndersonDarlingUnbinned = true;

const bool isOnChi2                  = true;
const bool isOnCramerVonMisesBinned  = false; 
const bool isOnAndersonDarlingBinned = false; 

const bool isCombinedAcceptanceSuiteTest = true;


void calculatePower() {

  std::cout << " Start of  calculatePower() " << std::endl
            << "\t isOnKolmogorovSmirnov = " << isOnKolmogorovSmirnov << std::endl
            << "\t isOnCramerVonMisesUnbinned = " << isOnCramerVonMisesUnbinned 
	    << std::endl
            << "\t isOnAndersonDarlingUnbinned = " << isOnAndersonDarlingUnbinned
            << std::endl
            << "\t isOnChi2 = " << isOnChi2 << std::endl
            << "\t isOnCramerVonMisesBinned = " << isOnCramerVonMisesBinned << std::endl
            << "\t isOnAndersonDarlingBinned = " << isOnAndersonDarlingBinned 
	    << std::endl
            << "\t isCombinedAcceptanceSuiteTest = " << isCombinedAcceptanceSuiteTest
	    << std::endl;

  // --- AIDA initialization ---

  // Creating the analysis factory.
  std::auto_ptr<AIDA::IAnalysisFactory> af( AIDA_createAnalysisFactory() );

  // Creating the tree factory.
  std::auto_ptr<AIDA::ITreeFactory> tf( af->createTreeFactory() );

  // Creating a tree in memory.
  std::auto_ptr<AIDA::ITree> tree( tf->create() );

  // Creating a histogram factory, whose histograms will be handled by the tree.
  std::auto_ptr<AIDA::IHistogramFactory> hf( af->createHistogramFactory( *tree ) );
  
  // Creating the two 1D clouds once, and then reset/fill
  // them for each pair of distributions which must be compared.
  AIDA::ICloud1D& firstDistribution  = *( hf->createCloud1D( "firstCloud" ) );
  AIDA::ICloud1D& secondDistribution = *( hf->createCloud1D( "secondCloud" ) );

  // --- Now the real code test ---

  const int numberPossibleStatisticalTests = 7; // 6 + 1 combined.
  StatisticsComparator< KolmogorovSmirnovComparisonAlgorithm >       comparatorKS; 
  StatisticsComparator< Chi2ComparisonAlgorithm >                    comparatorC2;
  StatisticsComparator< CramerVonMisesUnbinnedComparisonAlgorithm >  comparatorCM;
  StatisticsComparator< CramerVonMisesBinnedComparisonAlgorithm >    comparatorCMB;
  StatisticsComparator< AndersonDarlingUnbinnedComparisonAlgorithm > comparatorAD;
  StatisticsComparator< AndersonDarlingBinnedComparisonAlgorithm >   comparatorADB;

  std::vector< StatisticsComparatorStruct > comparatorVec;
  bool isAtLeastOneBinnedComparatorOn = false;
  bool isAtLeastOneUnbinnedComparatorOn = false;

  for ( int iComparatorCase = 0; 
	iComparatorCase < numberPossibleStatisticalTests; iComparatorCase++ ){
  
    bool isOnTest = false;
    bool isBinnedTest;
    std::string nameComparator;
    switch ( iComparatorCase ) {
    case 0 : {
      if ( isOnKolmogorovSmirnov ) {
	isOnTest = true;
        isBinnedTest = false;
	isAtLeastOneUnbinnedComparatorOn = true;
	nameComparator = "\t *** Kolmogorov-Smirnov        test *** ";   
      }
      break;
    }
    case 1 : {
      if ( isOnChi2 ) {
	isOnTest = true;
        isBinnedTest = true;
	isAtLeastOneBinnedComparatorOn = true;
	nameComparator = "\t *** Chi2                      test *** "; 
      }
      break;
    }
    case 2 : {
      if ( isOnCramerVonMisesUnbinned ) {
	isOnTest = true;
        isBinnedTest = false;
	isAtLeastOneUnbinnedComparatorOn = true;
	nameComparator = "\t *** Cramer-Von Mises Unbinned test *** "; 
      }
      break;
    }
    case 3 : {
      if ( isOnCramerVonMisesBinned ) {
	isOnTest = true;
        isBinnedTest = true;
	isAtLeastOneBinnedComparatorOn = true;
	nameComparator = "\t *** Cramer-Von Mises Binned   test *** "; 
      }
      break;
    }
    case 4 : {
      if ( isOnAndersonDarlingUnbinned ) {
	isOnTest = true;
        isBinnedTest = false;
	isAtLeastOneUnbinnedComparatorOn = true;
	nameComparator = "\t *** Anderson-Darling Unbinned test *** "; 
      }
      break;
    }
    case 5 : {
      if ( isOnAndersonDarlingBinned ) {
	isOnTest = true;
        isBinnedTest = true;
	isAtLeastOneBinnedComparatorOn = true;
	nameComparator = "\t *** Anderson-Darling Binned   test *** "; 
      }
      break;
    }
    case 6 : {
      if ( isCombinedAcceptanceSuiteTest ) {
	isOnTest = true;
        isBinnedTest = true;
	isAtLeastOneBinnedComparatorOn = true;
	nameComparator = "\t *** Combined Acceptance Suite test *** "; 
      }
      break;
    }
    default : std::cout << " ***WRONG*** it should never happen! " << std::endl;     
    }
    if ( isOnTest ) {
      StatisticsComparatorStruct theCompStruct;
      theCompStruct.id = iComparatorCase; 
      theCompStruct.fullName = nameComparator;
      theCompStruct.isBinned = isBinnedTest;
      comparatorVec.push_back( theCompStruct );    
    }

  }

  if ( ! comparatorVec.size() ) {
    std::cout << " === NONE Statistics Test IS ACTIVE === " << std::endl;
    return;
  }

  PowerCalculator powerCalculator;
  powerCalculator.setIsBinnedTest( isAtLeastOneBinnedComparatorOn );

  int countTotalNumberOfTests = 0;
  int countNumberOfTestsWonByKS = 0; // Smallest p-values won by Kolmogorov-Smirnov 
  int countNumberOfTestsWonByC2 = 0; //    "         "     "   " Chi2
  int countNumberOfTestsWonByCM = 0; //    "         "     "   " Cramer-von Mises
  int countNumberOfTestsWonByAD = 0; //    "         "     "   " Anderson-Darling

  // Loop over the pseudoexperiments.
  for ( int iPseudoExp = 0; 
	iPseudoExp < PowerCalculator::numberPseudoExps; iPseudoExp++ ) {

    // std::cout << " Starting Loop PseudoExp = " << iPseudoExp 
    //           << std::endl; //***DEBUG***
 
    std::vector< GeneratedPseudoExpStruct > generation; 
    powerCalculator.generateOnePseudoExp( generation );
    // Generate the  iComparsion-th  pseudoexperiment.

    std::vector< TestResultPseudoExpStruct > resultVec;
	
    // Loop over the two samples combinations that have been considered.
    for ( std::vector< GeneratedPseudoExpStruct >::iterator itGen = generation.begin();
          itGen != generation.end(); ++itGen ) {
 
      // std::cout << "\t Starting Loop Sample Combinations " << std::endl; //***DEBUG***

      std::vector< double > & firstVector  = itGen->firstSample;
      std::vector< double > & secondVector = itGen->secondSample;
      std::vector< double > & binsVector   = itGen->binning;

      AIDA::IHistogram1D *ptrHisto1 = 0, *ptrHisto2 = 0; 
	
      if ( isAtLeastOneBinnedComparatorOn ) {

	if ( binsVector.size() <= 1 ) {
	  std::cout << "***WRONG*** No bins! " << std::endl;
	  return;
	}
	// Creating a 1D Histogram
	ptrHisto1 = hf->createHistogram1D( "firstHisto", "The first histogram", 
					   binsVector );
	ptrHisto2 = hf->createHistogram1D( "secondHisto", "The second histogram", 
					   binsVector );

	// Fill the histograms with the vectors.
	for ( std::vector< double >::iterator itVec1 = firstVector.begin();
	      itVec1 != firstVector.end(); ++itVec1 ) {
	  ptrHisto1->fill( *itVec1 );
	}
	for ( std::vector< double >::iterator itVec2 = secondVector.begin();
	      itVec2 != secondVector.end(); ++itVec2 ) {
	  ptrHisto2->fill( *itVec2 );
	}

	//***DEBUG***
	// std::cout << "\t \t \t firstHistogram \t secondHistogram" << std::endl
	// 	  << "\t \t isFixedBinning = " 
	// 	  << ptrHisto1->axis().isFixedBinning() << std::endl
	// 	  << "\t \t bins : " << ptrHisto1->axis().bins() 
	// 	  << "\t edges : " << ptrHisto1->axis().lowerEdge()
	// 	  << "\t" << ptrHisto1->axis().upperEdge() << std::endl
	//           << "\t \t allEntries = " 
	//           << ptrHisto1->allEntries() << "\t \t"
	//           << ptrHisto2->allEntries() << std::endl
	//           << "\t \t inside entries = " 
	//           << ptrHisto1->entries() << "\t \t"
	//           << ptrHisto2->entries() << std::endl
	//           << "\t \t underflow+overflow entries = " 
	//           << ptrHisto1->extraEntries() << "\t \t"
	//           << ptrHisto2->extraEntries() << std::endl
	//           << "\t \t mean = " 
	//           << ptrHisto1->mean() << "\t \t"
	//           << ptrHisto2->mean() << std::endl
	//           << "\t \t rms = " 
	//           << ptrHisto1->rms() << "\t \t"
	//           << ptrHisto2->rms() << std::endl;
	// for ( int ii = 0 ; ii < ptrHisto1->axis().bins(); ii++ ) {
	//   std::cout << "\t \t bin = " << ii 
	// 	    << "\t edges : " << ptrHisto1->axis().binLowerEdge( ii )
	// 	    << "\t \t " << ptrHisto1->axis().binUpperEdge( ii ) << std::endl
	//             << "\t \t \t entries = " 
	// 	    << ptrHisto1->binEntries( ii ) << "\t \t"
	// 	    << ptrHisto2->binEntries( ii ) << std::endl
	// 	    << "\t \t \t weights = " 
	// 	    << ptrHisto1->binHeight( ii ) << "\t \t"
	// 	    << ptrHisto2->binHeight( ii ) << std::endl
	// 	    << "\t \t \t errors = " 
	// 	    << ptrHisto1->binError( ii ) << "\t \t"
	// 	    << ptrHisto2->binError( ii ) << std::endl;
	// }
	//***endDEBUG***

      }

      if ( isAtLeastOneUnbinnedComparatorOn ) {
	  
	// Reset the two 1D clouds, and fill them with the vectors.
	firstDistribution.reset();
	for ( std::vector< double >::iterator jtVec1 = firstVector.begin();
	      jtVec1 != firstVector.end(); ++jtVec1 ) {
	  firstDistribution.fill( *jtVec1 );
	}
	secondDistribution.reset();
	for ( std::vector< double >::iterator jtVec2 = secondVector.begin();
	      jtVec2 != secondVector.end(); ++jtVec2 ) {
	  secondDistribution.fill( *jtVec2 );
	}
	
	//***DEBUG***
	// std::cout << "\t firstDistribution  entries=" 
	//           << firstDistribution.entries() << std::endl; 
	// for ( int i = 0 ; i < firstDistribution.entries(); i++ ) {
	//   std::cout << "\t \t" << i << "\t" << firstDistribution.value( i ) 
	//             << std::endl;
	// }
	// std::cout << "\t secondDistribution  entries=" 
	//           << secondDistribution.entries() << std::endl; 
	// for ( int j = 0 ; j < secondDistribution.entries(); j++ ) {
	//   std::cout << "\t \t" << j << "\t" << secondDistribution.value( j ) 
	//             << std::endl;
	// }
	//***endDEBUG***
	
      }

      // Loop over the statistics tests.
      for ( std::vector< StatisticsComparatorStruct >::iterator jtStat = 
	      comparatorVec.begin(); jtStat != comparatorVec.end(); ++jtStat ) { 

	// std::cout << "\t \t Starting Loop Statistics Tests " 
        //           << std::endl; //***DEBUG***

	StatisticsTesting::ComparisonResult result;

	// Statistical test. 
	switch ( jtStat->id ) {
	case 0 : {
	  result = comparatorKS.compare( firstDistribution, secondDistribution );
	  break;
	}
	case 1 : {
	  result = comparatorC2.compare( *ptrHisto1, *ptrHisto2 );
	  break;
	}
	case 2 : {
	  result = comparatorCM.compare( firstDistribution, secondDistribution );
	  break;
	}
	case 3 : {
	  result = comparatorCMB.compare( *ptrHisto1, *ptrHisto2 );
	  break;
	}
	case 4 : {
	  result = comparatorAD.compare( firstDistribution, secondDistribution );
	  break;
	}
	case 5 : {
	  result = comparatorADB.compare( *ptrHisto1, *ptrHisto2 );
	  break;
	}
	case 6 : {
	  // The result is the one corresponding to the statistical test
	  // with the lowest p-value between the tests used in the 
          // Acceptance Suite: Chi2, Kolmogorov-Smirnov, Cramer-von Mises, 
          // Anderson-Darling.
	  StatisticsTesting::ComparisonResult resultKS = 
	    comparatorKS.compare( firstDistribution, secondDistribution );
	  StatisticsTesting::ComparisonResult resultC2 = 
	    comparatorC2.compare( *ptrHisto1, *ptrHisto2 );
	  StatisticsTesting::ComparisonResult resultCM = 
	    comparatorCM.compare( firstDistribution, secondDistribution );
	  StatisticsTesting::ComparisonResult resultAD = 
	    comparatorAD.compare( firstDistribution, secondDistribution );
	  int winningCase = 0;
          result = resultKS;
	  if ( resultC2.quality() < result.quality() ) {
	    result = resultC2;
	    winningCase = 1;
	  }
	  if ( resultCM.quality() < result.quality() ) {
	    result = resultCM;
	    winningCase = 2;
	  }
	  if ( resultAD.quality() < result.quality() ) {
	    result = resultAD;
	    winningCase = 4;
	  }
	  switch ( winningCase ) {
	    case 0 : countNumberOfTestsWonByKS++; break;
	    case 1 : countNumberOfTestsWonByC2++; break;
	    case 2 : countNumberOfTestsWonByCM++; break;
	    case 4 : countNumberOfTestsWonByAD++; break;
          }
	  countTotalNumberOfTests++;
	  break;
	}
	default : std::cout << " ***WRONG*** it should never happen! " << std::endl;     
	}
      
	struct TestResultPseudoExpStruct theResultStruct;
        theResultStruct.pseudoExpNum = iPseudoExp;
        theResultStruct.idParents    = itGen->idParents;
        theResultStruct.idStatTest   = jtStat->id;
        theResultStruct.distance     = result.distance();
        theResultStruct.pValue       = result.quality();
        theResultStruct.ndf          = result.ndf();

        resultVec.push_back( theResultStruct );

	//std::cout << "\t PseudoExp=" << theResultStruct.pseudoExpNum
	//	    << "  parents="    << theResultStruct.idParents
	// 	    << "  statTest="   << theResultStruct.idStatTest
	// 	    << "  d="          << theResultStruct.distance 
	// 	    << "  p="          << theResultStruct.pValue 
	// 	    << "  ndf="        << theResultStruct.ndf
	// 	    << std::endl;        //***DEBUG***
	
      } // End loop over the statistics tests.

      // Remove the histograms from the tree.
      tree->rm( "firstHisto" );
      tree->rm( "secondHisto" );

    } // End loop over the two samples combinations.

    powerCalculator.storeTestResultPseudoExp( resultVec );

  } // End loop over pseudoexperiments.
    
  powerCalculator.printResults( comparatorVec );

  if ( countTotalNumberOfTests > 0 ) {
    std::cout << " --- Winning tests for the smallest p-values --- " 
	      << std::endl
	      << " KS : " 
	      << 100.0 * ( static_cast<double>( countNumberOfTestsWonByKS ) / 
			   static_cast<double>( countTotalNumberOfTests ) ) 
	      << " %" << std::endl
	      << " C2 : " 
	      << 100.0 * ( static_cast<double>( countNumberOfTestsWonByC2 ) / 
			   static_cast<double>( countTotalNumberOfTests ) ) 
	      << " %" << std::endl
	      << " CM : " 
	      << 100.0 * ( static_cast<double>( countNumberOfTestsWonByCM ) / 
			   static_cast<double>( countTotalNumberOfTests ) ) 
	      << " %" << std::endl
	      << " AD : " 
	      << 100.0 * ( static_cast<double>( countNumberOfTestsWonByAD ) / 
			   static_cast<double>( countTotalNumberOfTests ) ) 
	      << " %" << std::endl
              << " ----------------------------------------------- "
              << std::endl;
  }

}


//==================================================================================
// MAIN
//==================================================================================

int main (int, char **) {

  //***LOOKHERE*** Change the random seed here if you want.
  RanluxEngine defaultEngine( 1234, 4 ); 
  // RanluxEngine defaultEngine( 56789, 4 ); 
  
  HepRandom::setTheEngine( &defaultEngine );
  
  std::cout << std::endl << std::endl
	    << " ================ power : STARTS ================= "
	    << std::endl << std::endl;
  
  calculatePower();
  
  std::cout << std::endl << std::endl
	    << " ================ power : ENDS ================= "
	    << std::endl << std::endl;
  
}


//==================================================================================
//==================================================================================
