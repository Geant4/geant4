//========================================================================
//                         Program  pvalue.cpp
//                         -------------------
// This program reads in input two HBOOK files: 
//          1)  ntuple_a.hbook 
//          2)  ntuple_b.hbook
// which contain each a ntuple, containing physics information about
// the showering of a certain incoming beam particle into a sampling 
// calorimeter. One of the ntuple can be obtained with a certain version
// of Geant4, for instance 6.0, and the other with another version,
// for instance 6.1 . 
// The goal of this program is to compare the two ntuples, in such a 
// way to do "regression testing", i.e. to see if the two versions are
// statistically compatible, or there is at least one observable which
// is significantly different, and that needs to be understood (it can
// be a bug, or an improvement, or a new feature). To do so, the program
// projects the ntuples into different histograms, for binned statistical
// tests, and into "clouds" (AIDA objects which are basically vectors of
// values, i.e. a kind of "unbinned histogram"), for unbinned statistical
// tests. The output of this program consists of:
//          a)  cloudsA.xml, cloudsB.xml : XML files containing the 
//                                         clouds obtained from the
//              first ntuple (ntuple_a.hbook) and the second one (ntuple_b.hbook),
//              respectively.
//          b)  histo_a.hbook, histo_b.hbook : HBOOK files containing
//                                             the histograms obtained
//              from the first ntuple (ntuple_a.hbook) and the second one
//              (ntuple_b.hbook), respectively.
//          c)  printing on the screen the results of the statistical
//              tests, i.e. the probability, called "p-value", that 
//              two histograms (in the case of binned statistics test),
//              or two clouds (in the case of unbinned statistics test)
//              are coming from the same parent distribution (whatever
//              it is). 
//
//========================================================================

#include <memory>
#include <iostream>
#include <iomanip>

#include "AIDA/AIDA.h"
#include "AIDA/ICloud1D.h" 
#include "StatisticsTesting/StatisticsComparator.h"
#include "ComparisonResult.h"

#include "Chi2ComparisonAlgorithm.h"                    
//#include "CramerVonMisesBinnedComparisonAlgorithm.h"  
//#include "AndersonDarlingBinnedComparisonAlgorithm.h" 

#include "KolmogorovSmirnovComparisonAlgorithm.h"         
// #include "CramerVonMisesUnbinnedComparisonAlgorithm.h" 
// #include "AndersonDarlingUnbinnedComparisonAlgorithm.h"

#include <vector>
#include <sstream>

using namespace StatisticsTesting; 

const double pvalueThreshold = 0.01;   // A p-value below this threshold is
                                       // considered as a "warning" that the
                                       // two distributions are not coming
                                       // from the same parent.

//==================
// GLOBAL FUNCTIONS
//==================

// This function converts an integer into a string.
std::string toString( int n ) {
  std::ostringstream buf;
  buf << n << std::ends;      // ends adds the end of string character \0 .
  return buf.str();
}


// The following function has in input an histogram, and it fills
// a vector with the bin edges of that histogram.
void computeBinEdges( const AIDA::IHistogram1D & histogramIn ,   // input
                      std::vector< double > & binEdgesOut ) {    // output
  const AIDA::IAxis& xAxis = histogramIn.axis();
  for ( int iBin = 0; iBin < xAxis.bins(); ++iBin ) {
    binEdgesOut.push_back( xAxis.binLowerEdge( iBin ) );
  }
  binEdgesOut.push_back( xAxis.binUpperEdge( xAxis.bins() - 1 ) ); 
  return;
}


// The following function has in input two clouds, and it converts
// them to histograms, with the same binning. The histograms are
// kept inside the clouds themselves. The function returns true
// (false) is the conversion succeed (fail).
bool convertCloudsToTheSameHistogram( AIDA::ICloud1D & cloud1, 
				      AIDA::ICloud1D & cloud2 ) {
  bool succeed = true;
  AIDA::ICloud1D & cloudA = cloud1;
  AIDA::ICloud1D & cloudB = cloud2;

  // Decide which is the cloud which should follow the binning
  // of the other one: cloudA is leading, cloudB is following.
  // The algorithm is simple: cloudA is always the cloud1 
  // unless cloud1 is wider on both sides than cloud2.
  if ( cloud1.lowerEdge() < cloud2.lowerEdge()  &&
       cloud1.upperEdge() > cloud2.upperEdge() ) {
    cloudA = cloud2;
    cloudB = cloud1;
  }    

  // Force the conversion of cloudB to the same binning as the 
  // cloudA's natural conversion to an histogram. 
  if ( ! cloudA.convertToHistogram() ) {
    succeed = false;
  } else {
    std::vector<double> binEdges;
    computeBinEdges( cloudA.histogram(), binEdges );
    if ( ! cloudB.convert( binEdges ) ) {
      succeed = false;
    }
  }    

  if ( ! succeed ) {
    std::cerr << "***ERROR*** : Cannot convert a cloud into histogram " << std::endl;
  }
  return succeed;
}


//======
// MAIN
//======

int main (int, char **) {

  std::cout << " START " << std::endl;

  // ----------------------------------------------------------
  // AIDA initialization, and reading of the two input ntuples.
  // ----------------------------------------------------------

  // Creating the analysis factory.
  std::auto_ptr<AIDA::IAnalysisFactory> af( AIDA_createAnalysisFactory() );

  // Creating the tree factory.
  std::auto_ptr<AIDA::ITreeFactory> tf( af->createTreeFactory() );
 
  // Read the two trees from the two input files.
  bool readOnly  = true;  
  bool createNew = false; 
  std::auto_ptr<AIDA::ITree> treeHBookA( tf->create( "ntuple_a.hbook", "hbook", 
						     readOnly, createNew ) );
  if ( ! treeHBookA.get() ) {
    std::cout << " ERROR: unable to open file  ntuple_a.hbook " << std::endl;
    return -1;
  } else {
    // std:: cout << " OK : opened file  ntuple_a.hbook " << std::endl   //***DEBUG***
    //           << "\t The content is: " << treeHBookA->ls() << std::endl;
  }

  std::auto_ptr<AIDA::ITree> treeHBookB( tf->create( "ntuple_b.hbook", "hbook", 
                                                     readOnly, createNew ) );
  if ( ! treeHBookB.get() ) {
    std::cout << " ERROR: unable to open file  ntuple_b.hbook " << std::endl;
    return -1;
  } else {
    // std:: cout << " OK : opened file  ntuple_b.hbook " << std::endl   //***DEBUG***
    //            << "\t The content is: " << treeHBookB->ls() << std::endl;

  }

  // std::cout << " treeHBookA " << (int) treeHBookA.get() << std::endl; //***DEBUG***
  AIDA::ITuple* p_tpA = dynamic_cast<AIDA::ITuple*>(  ( treeHBookA->find( "/1" ) ) );
  if ( ! p_tpA ) { 
     std::cerr << "Error finding /hbookA/1 " << std::endl; 
     return -1;
  }
  // std::cout << (int) p_tpA << std::endl; //***DEBUG*** 
    
  AIDA::ITuple& tpA = *p_tpA; 

  //***DEBUG***
  // std::cout << "Tuple A title : " << tpA.title() << std::endl;
  // std::cout << "Tuple A variables : " << std::endl;
  // for ( int i = 0; i < tpA.columns(); ++i ) {
  //   std::cout << tpA.columnName(i) << "\t" << tpA.columnType(i) << std::endl;
  // }

  // ---------------------------------------------------------------
  // Get the number of elements of the vector L and R of the ntuple.
  // ---------------------------------------------------------------
  int i_nLayers = tpA.findColumn( "nLayers" ); 
  int i_nBinR = tpA.findColumn( "nBinR" ); 
  tpA.start(); tpA.next();
  const int numberOfReplicas   = tpA.getInt( i_nLayers ); 
  const int numberOfRadiusBins = tpA.getInt( i_nBinR ); 

  std::cout << " numberOfReplicas   = " << numberOfReplicas << std::endl
            << " numberOfRadiusBins = " << numberOfRadiusBins 
	    << std::endl;                                       //***DEBUG***

  // -------------------------------------------------------------
  // Creating the clouds, and fill them by projecting the ntuples.
  // -------------------------------------------------------------

  // Creating two trees on two files where to save all clouds.
  readOnly  = false;  
  createNew = true; 
  std::auto_ptr<AIDA::ITree> treeCloudsA( tf->create("cloudsA.xml", "xml", 
						     readOnly, createNew) );
  std::auto_ptr<AIDA::ITree> treeCloudsB( tf->create("cloudsB.xml", "xml", 
						     readOnly, createNew) );

  // Creating two histogram factories attached to the two trees.
  std::auto_ptr<AIDA::IHistogramFactory> 
    hfA( af->createHistogramFactory( *treeCloudsA ) );
  std::auto_ptr<AIDA::IHistogramFactory> 
    hfB( af->createHistogramFactory( *treeCloudsB ) );

  // Create the clouds.
  AIDA::ICloud1D * cA1 = hfA->createCloud1D( "Energy deposit in Active layers, A" ); 
  AIDA::ICloud1D * cA2 = hfA->createCloud1D( "Energy deposit in All layers, A" ); 
  std::vector< AIDA::ICloud1D * > cAL;   
  std::vector< AIDA::ICloud1D * > cAR;   
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    std::string name = "LA" + toString( i );
    cAL.push_back( hfA->createCloud1D( name ) ); 
    // std::cout << "\t Created cloud: name = " << name << std::endl; //***DEBUG***
  }
  for ( int i = 0; i < numberOfRadiusBins; i++ ) {
    std::string name = "RA" + toString( i );
    cAR.push_back( hfA->createCloud1D( name ) ); 
    // std::cout << "\t Created cloud: name = " << name << std::endl;  //***DEBUG***
  }

  int iEdepAct = tpA.findColumn( "EDEP_ACT" ); 
  int iEdepCal = tpA.findColumn( "EDEP_CAL" ); 
  int iL = tpA.findColumn( "L" );
  int iR = tpA.findColumn( "R" ); 

  // Fill the clouds for the first ntuple.
  tpA.start(); 
  while ( tpA.next() ) { 
    float eDepAct = tpA.getFloat( iEdepAct ); 
    cA1->fill( eDepAct );
    // std::cout << " Filling cloud A1 with " << eDepAct << std::endl; //***DEBUG***
    float eDepCal = tpA.getFloat( iEdepCal ); 
    cA2->fill( eDepCal ); 
    // std::cout << " Filling cloud A2 with " << eDepCal << std::endl; //***DEBUG***
    AIDA::ITuple& tL = *( tpA.getTuple( iL ) );
    tL.start(); 
    int indL = tL.findColumn( "L" );
    for ( int i = 0; i < tL.rows() ; i++ ) { 
      tL.setRow( i ); 
      float value = tL.getFloat( indL ); 
      cAL[ i ]->fill( value );
      // std::cout << " Filling cloud AL[" << i << "] with " << value 
      //           << std::endl; //***DEBUG***
    }  
    AIDA::ITuple& tR = *( tpA.getTuple( iR ) );
    tR.start(); 
    int indR = tR.findColumn( "R" );
    for ( int i = 0; i < tR.rows() ; i++ ) { 
      tR.setRow( i ); 
      float value = tR.getFloat( indR ); 
      cAR[ i ]->fill( value );
      // std::cout << " Filling cloud AR[" << i << "] with " << value 
      //           << std::endl; //***DEBUG***
    }  
  }

  AIDA::ITuple& tpB = dynamic_cast<AIDA::ITuple&>( * ( treeHBookB->find( "/1" ) ) );

  //***DEBUG***
  // std::cout << "Tuple B title : " << tpB.title() << std::endl;
  // std::cout << "Tuple B variables : " << std::endl;
  // for ( int i = 0; i < tpB.columns(); ++i ) {
  //   std::cout << tpB.columnName(i) << "\t" << tpB.columnType(i) << std::endl;
  // }

  AIDA::ICloud1D * cB1 = hfB->createCloud1D( "Energy deposit in Active layers, B" ); 
  AIDA::ICloud1D * cB2 = hfB->createCloud1D( "Energy deposit in All layers, B" ); 
  std::vector< AIDA::ICloud1D * > cBL;   
  std::vector< AIDA::ICloud1D * > cBR;   
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    std::string name = "LB" + toString( i );
    cBL.push_back( hfB->createCloud1D( name ) ); 
    // std::cout << "\t Created cloud: name = " << name << std::endl; //***DEBUG***
  }
  for ( int i = 0; i < numberOfRadiusBins; i++ ) {
    std::string name = "RB" + toString( i );
    cBR.push_back( hfB->createCloud1D( name ) ); 
    // std::cout << "\t Created cloud: name = " << name << std::endl;  //***DEBUG***
  }

  // Fill the clouds for the second ntuple.
  tpB.start(); 
  while ( tpB.next() ) { 
    float eDepAct = tpB.getFloat( iEdepAct ); 
    cB1->fill( eDepAct );
    // std::cout << " Filling cloud B1 with " << eDepAct << std::endl; //***DEBUG***
    float eDepCal = tpB.getFloat( iEdepCal ); 
    cB2->fill( eDepCal );
    // std::cout << " Filling cloud B2 with " << eDepCal << std::endl; //***DEBUG***
    AIDA::ITuple& tL = *( tpB.getTuple( iL ) );
    tL.start(); 
    int indL = tL.findColumn( "L" );
    for ( int i = 0; i < tL.rows() ; i++ ) { 
      tL.setRow( i ); 
      float value = tL.getFloat( indL ); 
      cBL[ i ]->fill( value );
      // std::cout << " Filling cloud BL[" << i << "] with " << value 
      //           << std::endl; //***DEBUG***
    }  
    AIDA::ITuple& tR = *( tpB.getTuple( iR ) );
    tR.start(); 
    int indR = tR.findColumn( "R" );
    for ( int i = 0; i < tR.rows() ; i++ ) { 
      tR.setRow( i ); 
      float value = tR.getFloat( indR ); 
      cBR[ i ]->fill( value );
      // std::cout << " Filling cloud BR[" << i << "] with " << value 
      //           << std::endl; //***DEBUG***
    }  
  }

  //***DEBUG*** : Check the content of the clouds.
  // std::cout << " Debugging info about cloud 1 :  A \t B " << std::endl
  //           << "\t sumOfWeights : " << cA1->sumOfWeights() 
  // 	    << "\t"                 << cB1->sumOfWeights() << std::endl
  //           << "\t lowerEdge    : " << cA1->lowerEdge() 
  // 	    << "\t"                 << cB1->lowerEdge() << std::endl
  //           << "\t upperEdge    : " << cA1->upperEdge() 
  // 	    << "\t"                 << cB1->upperEdge() << std::endl
  //           << "\t mean         : " << cA1->mean() 
  // 	    << "\t"                 << cB1->mean() << std::endl
  //           << "\t rms          : " << cA1->rms() 
  // 	    << "\t"                 << cB1->rms() << std::endl;
  // for ( int i = 0; i < cA1->sumOfWeights(); i++ ) {
  //   std::cout << "\t i=" << i << "\t" << cA1->value( i ) 
  //             << "\t"                 << cB1->value( i ) << std::endl;
  // } 
  // std::cout << " Debugging info about cloud 2 :  A \t B " << std::endl
  //           << "\t sumOfWeights : " << cA2->sumOfWeights() 
  // 	    << "\t"                 << cB2->sumOfWeights() << std::endl
  //           << "\t lowerEdge    : " << cA2->lowerEdge() 
  // 	    << "\t"                 << cB2->lowerEdge() << std::endl
  //           << "\t upperEdge    : " << cA2->upperEdge() 
  // 	    << "\t"                 << cB2->upperEdge() << std::endl
  //           << "\t mean         : " << cA2->mean() 
  // 	    << "\t"                 << cB2->mean() << std::endl
  //           << "\t rms          : " << cA2->rms() 
  // 	    << "\t"                 << cB2->rms() << std::endl;
  // for ( int i = 0; i < numberOfReplicas; i++ ) {
  //   std::cout << " Debugging info about cloud L[" << i << "] :  A \t B " << std::endl
  // 	      << "\t sumOfWeights : " << cAL[i]->sumOfWeights() 
  // 	      << "\t"                 << cBL[i]->sumOfWeights() << std::endl
  // 	      << "\t lowerEdge    : " << cAL[i]->lowerEdge() 
  // 	      << "\t"                 << cBL[i]->lowerEdge() << std::endl
  // 	      << "\t upperEdge    : " << cAL[i]->upperEdge() 
  // 	      << "\t"                 << cBL[i]->upperEdge() << std::endl
  // 	      << "\t mean         : " << cAL[i]->mean() 
  // 	      << "\t"                 << cBL[i]->mean() << std::endl
  // 	      << "\t rms          : " << cAL[i]->rms() 
  // 	      << "\t"                 << cBL[i]->rms() << std::endl;
  // }
  // for ( int i = 0; i < numberOfRadiusBins; i++ ) {
  //   std::cout << " Debugging info about cloud R[" << i << "] :  A \t B " << std::endl
  // 	      << "\t sumOfWeights : " << cAR[i]->sumOfWeights() 
  // 	      << "\t"                 << cBR[i]->sumOfWeights() << std::endl
  // 	      << "\t lowerEdge    : " << cAR[i]->lowerEdge() 
  // 	      << "\t"                 << cBR[i]->lowerEdge() << std::endl
  // 	      << "\t upperEdge    : " << cAR[i]->upperEdge() 
  // 	      << "\t"                 << cBR[i]->upperEdge() << std::endl
  // 	      << "\t mean         : " << cAR[i]->mean() 
  // 	      << "\t"                 << cBR[i]->mean() << std::endl
  // 	      << "\t rms          : " << cAR[i]->rms() 
  // 	      << "\t"                 << cBR[i]->rms() << std::endl;
  // }
 
  // Flushing the clouds into the files
  treeCloudsA->commit();
  treeCloudsB->commit();
  
  // ------------------------------------------------------------
  // Creating the histograms, by converting copies of the clouds.
  // ------------------------------------------------------------

  // Copy the clouds (before converting them to histograms)
  // Notice that it is important that the names of the copies are
  // different than the original, because they belong to the same
  // histogram factories, and hence trees, and therefore the names
  // are used as identifiers.
  AIDA::ICloud1D * cA1bis = 
    hfA->createCopy( "Energy deposit in Active layers, A bis", *cA1 );
  AIDA::ICloud1D * cB1bis = 
    hfB->createCopy( "Energy deposit in Active layers, B bis", *cB1 );
  AIDA::ICloud1D * cA2bis = 
    hfA->createCopy( "Energy deposit in All layers, A bis", *cA2 ); 
  AIDA::ICloud1D * cB2bis = 
    hfB->createCopy( "Energy deposit in All layers, B bis", *cB2 ); 
  std::vector< AIDA::ICloud1D * > cALbis;
  std::vector< AIDA::ICloud1D * > cBLbis;
  std::vector< AIDA::ICloud1D * > cARbis;
  std::vector< AIDA::ICloud1D * > cBRbis;
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    std::string name = "LA" + toString( i ) + "bis";
    cALbis.push_back( hfA->createCopy( name, *cAL[i] ) ); 
    // std::cout << "\t Created cloud: name = " << name << std::endl;  //***DEBUG*** 
    name = "LB" + toString( i ) + "bis";
    cBLbis.push_back( hfB->createCopy( name, *cBL[i] ) ); 
    // std::cout << "\t Created cloud: name = " << name << std::endl;  //***DEBUG***
  }
  for ( int i = 0; i < numberOfRadiusBins; i++ ) {
    std::string name = "RA" + toString( i ) + "bis";
    cARbis.push_back( hfA->createCopy( name, *cAR[i] ) ); 
    // std::cout << "\t Created cloud: name = " << name << std::endl;  //***DEBUG*** 
    name = "RB" + toString( i ) + "bis";
    cBRbis.push_back( hfB->createCopy( name, *cBR[i] ) ); 
    // std::cout << "\t Created cloud: name = " << name << std::endl;  //***DEBUG*** 
  }

  // Convert the copies of the clouds into histograms:
  // notice that the pairs of histograms that must be compared
  // must have the same binning.

  if ( ! convertCloudsToTheSameHistogram( *cA1bis, *cB1bis )  ||
       ! convertCloudsToTheSameHistogram( *cA2bis, *cB2bis ) ) return -1;
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    if ( ! convertCloudsToTheSameHistogram( *cALbis[i], *cBLbis[i] ) ) return -1;
  }
  for ( int i = 0; i < numberOfRadiusBins; i++ ) {
    if ( ! convertCloudsToTheSameHistogram( *cARbis[i] , *cBRbis[i] ) ) return -1;
  }

  // Creating a tree mapped to an Hbook file where to store all the
  // histograms converted from clouds.
  readOnly = false;
  createNew = true;
  std::auto_ptr<AIDA::ITree> treeHBookAout( tf->create("histo_a.hbook", "hbook", 
						       readOnly, createNew ) );
  std::auto_ptr<AIDA::ITree> treeHBookBout( tf->create("histo_b.hbook", "hbook", 
						       readOnly, createNew ) );   
  // Creating a histogram factory, whose histograms will be handled by the tree
  std::auto_ptr<AIDA::IHistogramFactory> 
    hfacA( af->createHistogramFactory( *treeHBookAout ) );
  std::auto_ptr<AIDA::IHistogramFactory> 
    hfacB( af->createHistogramFactory( *treeHBookBout ) );

  hfacA->createCopy( "11", cA1bis->histogram() );
  hfacB->createCopy( "11", cB1bis->histogram() );
  hfacA->createCopy( "12", cA2bis->histogram() );
  hfacB->createCopy( "12", cB2bis->histogram() );
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    int idNum = i + 100;
    std::string idString = toString( idNum );
    hfacA->createCopy( idString, cALbis[i]->histogram() );
    hfacB->createCopy( idString, cBLbis[i]->histogram() );
  }
  for ( int i = 0; i < numberOfRadiusBins; i++ ) {
    int idNum = i + 200;
    std::string idString = toString( idNum );
    hfacA->createCopy( idString, cARbis[i]->histogram() );
    hfacB->createCopy( idString, cBRbis[i]->histogram() );
  }
  
  //***DEBUG*** : Check the content of the histograms.
  // std::cout << " Debugging info about histogram 1 :  A \t B " << std::endl
  //           << "\t allEntries   : " << cA1bis->histogram().allEntries() 
  // 	    << "\t"                 << cB1bis->histogram().allEntries() << std::endl
  //           << "\t extraEntries : " << cA1bis->histogram().extraEntries() 
  // 	    << "\t"                 << cB1bis->histogram().extraEntries() << std::endl
  //           << "\t mean         : " << cA1bis->histogram().mean() 
  // 	    << "\t"                 << cB1bis->histogram().mean() << std::endl
  //           << "\t rms          : " << cA1bis->histogram().rms() 
  // 	    << "\t"                 << cB1bis->histogram().rms() << std::endl;
  // Printing the contents of the histogram
  // std::cout << "\t X value \t entries \t \t Y value (height)" << std::endl;
  // const AIDA::IAxis& xAxis = cA1bis->histogram().axis();
  // for ( int iBin = 0; iBin < xAxis.bins(); ++iBin ) {
  //   std::cout << "\t"    << cA1bis->histogram().binMean( iBin )
  //             << "\t \t" << cA1bis->histogram().binEntries( iBin )
  //             << "\t"    << cB1bis->histogram().binEntries( iBin )
  //             << "\t \t" << cA1bis->histogram().binHeight( iBin )
  //             << "\t"    << cB1bis->histogram().binHeight( iBin )
  //             << std::endl;
  // }
  // std::cout << " Debugging info about histogram 2 :  A \t B " << std::endl
  //           << "\t allEntries   : " << cA2bis->histogram().allEntries() 
  // 	    << "\t"                 << cB2bis->histogram().allEntries() << std::endl
  //           << "\t extraEntries : " << cA2bis->histogram().extraEntries() 
  // 	    << "\t"                 << cB2bis->histogram().extraEntries() << std::endl
  //           << "\t mean         : " << cA2bis->histogram().mean() 
  // 	    << "\t"                 << cB2bis->histogram().mean() << std::endl
  //           << "\t rms          : " << cA2bis->histogram().rms() 
  // 	    << "\t"                 << cB2bis->histogram().rms() << std::endl;
  // for ( int i = 0; i < numberOfReplicas; i++ ) {
  //   std::cout << " Debugging info about histogram L[" << i << "] :  A \t B " << std::endl
  // 	      << "\t allEntries   : " << cALbis[i]->histogram().allEntries() 
  // 	      << "\t"                 << cBLbis[i]->histogram().allEntries() << std::endl
  // 	      << "\t extraEntries : " << cALbis[i]->histogram().extraEntries() 
  // 	      << "\t"                 << cBLbis[i]->histogram().extraEntries() << std::endl
  // 	      << "\t mean         : " << cALbis[i]->histogram().mean() 
  // 	      << "\t"                 << cBLbis[i]->histogram().mean() << std::endl
  // 	      << "\t rms          : " << cALbis[i]->histogram().rms() 
  // 	      << "\t"                 << cBLbis[i]->histogram().rms() << std::endl;
  // }
  // for ( int i = 0; i < numberOfRadiusBins; i++ ) {
  //   std::cout << " Debugging info about histogram R[" << i << "] :  A \t B " << std::endl
  // 	      << "\t allEntries   : " << cARbis[i]->histogram().allEntries() 
  // 	      << "\t"                 << cBRbis[i]->histogram().allEntries() << std::endl
  // 	      << "\t extraEntries : " << cARbis[i]->histogram().extraEntries() 
  // 	      << "\t"                 << cBRbis[i]->histogram().extraEntries() << std::endl
  // 	      << "\t mean         : " << cARbis[i]->histogram().mean() 
  // 	      << "\t"                 << cBRbis[i]->histogram().mean() << std::endl
  // 	      << "\t rms          : " << cARbis[i]->histogram().rms() 
  // 	      << "\t"                 << cBRbis[i]->histogram().rms() << std::endl;
  // }
  
  // Flushing the histograms into the file
  treeHBookAout->commit();
  treeHBookBout->commit();
  
  // ------------------------
  // Do the statistical tests
  // ------------------------

  StatisticsComparator< Chi2ComparisonAlgorithm > comparatorBin;                    
  // StatisticsComparator< CramerVonMisesBinnedComparisonAlgorithm > comparatorBin;
  // StatisticsComparator< AndersonDarlingBinnedComparisonAlgorithm > comparatorBin;

  StatisticsComparator< KolmogorovSmirnovComparisonAlgorithm > comparatorUnbin;
  // StatisticsComparator< CramerVonMisesUnbinnedComparisonAlgorithm > comparatorUnbin;
  // StatisticsComparator< AndersonDarlingUnbinnedComparisonAlgorithm > comparatorUnbin;
  
  std::cout << " ------------ Starting tests ------------- " << std::endl;
  
  std::cout << "Observable 1";
  ComparisonResult resultBin1 = comparatorBin.compare( cA1bis->histogram(), 
						       cB1bis->histogram() );
  ComparisonResult resultUnbin1 = comparatorUnbin.compare( *cA1, *cB1 ); 
  if ( resultBin1.quality() < pvalueThreshold  || 
       resultUnbin1.quality() < pvalueThreshold ) {
    std::cout << "\t ***WARNING***";
  }
  std::cout << std::endl;
  std::cout << "  Chi2" << "  d=" << resultBin1.distance()
	    << "  ndf=" << resultBin1.ndf() 
	    << "  pvalue=" << resultBin1.quality() << std::endl;
  std::cout << "  KS" << "  d=" << resultUnbin1.distance()
	    << "  ndf=" << resultUnbin1.ndf() 
	    << "  pvalue=" << resultUnbin1.quality() << std::endl;
      
  std::cout << "Observable 2";
  ComparisonResult resultBin2 = comparatorBin.compare( cA2bis->histogram(), 
						       cB2bis->histogram() );
  ComparisonResult resultUnbin2 = comparatorUnbin.compare( *cA2, *cB2 ); 
  if ( resultBin2.quality() < pvalueThreshold  || 
       resultUnbin2.quality() < pvalueThreshold ) {
    std::cout << "\t ***WARNING***";
  }
  std::cout << std::endl;
  std::cout << "  Chi2" << "  d=" << resultBin2.distance()
	    << "  ndf=" << resultBin2.ndf() 
	    << "  pvalue=" << resultBin2.quality() << std::endl;
  std::cout << "  KS" << "  d=" << resultUnbin2.distance()
            << "  ndf=" << resultUnbin2.ndf() 
            << "  pvalue=" << resultUnbin2.quality() << std::endl;
  
  std::vector< ComparisonResult > resultBinL;
  std::vector< ComparisonResult > resultUnbinL;
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    std::cout << "Observable L" << i;   
    ComparisonResult resultBinLvalue = comparatorBin.compare( cALbis[i]->histogram(),
  							      cBLbis[i]->histogram() );
    resultBinL.push_back( resultBinLvalue );
    ComparisonResult resultUnbinLvalue = comparatorUnbin.compare( *cAL[i], *cBL[i] );
    resultUnbinL.push_back( resultUnbinLvalue );
    if ( resultBinLvalue.quality() < pvalueThreshold ||
	 resultUnbinLvalue.quality() < pvalueThreshold ) {
      std::cout << "\t ***WARNING***";
    }
    std::cout << std::endl;
    std::cout << "  Chi2" << "  d=" << resultBinLvalue.distance()
  	      << "  ndf=" << resultBinLvalue.ndf() 
  	      << "  pvalue=" << resultBinLvalue.quality() << std::endl;
    std::cout << "  KS" << "  d=" << resultUnbinLvalue.distance()
  	      << "  ndf=" << resultUnbinLvalue.ndf() 
  	      << "  pvalue=" << resultUnbinLvalue.quality() << std::endl;
  }

  std::vector< ComparisonResult > resultBinR;
  std::vector< ComparisonResult > resultUnbinR;
  for ( int i = 0; i < numberOfRadiusBins; i++ ) {
    std::cout << "Observable R" << i;   
    ComparisonResult resultBinRvalue = comparatorBin.compare( cARbis[i]->histogram(),
  							      cBRbis[i]->histogram() );
    resultBinR.push_back( resultBinRvalue );
    ComparisonResult resultUnbinRvalue = comparatorUnbin.compare( *cAR[i], *cBR[i] );
    resultUnbinR.push_back( resultUnbinRvalue );
    if ( resultBinRvalue.quality() < pvalueThreshold ||
	 resultUnbinRvalue.quality() < pvalueThreshold ) {
      std::cout << "\t ***WARNING***";
    }
    std::cout << std::endl;
    std::cout << "  Chi2" << "  d=" << resultBinRvalue.distance()
  	      << "  ndf=" << resultBinRvalue.ndf() 
  	      << "  pvalue=" << resultBinRvalue.quality() << std::endl;
    std::cout << "  KS" << "  d=" << resultUnbinRvalue.distance()
  	      << "  ndf=" << resultUnbinRvalue.ndf() 
  	      << "  pvalue=" << resultUnbinRvalue.quality() << std::endl;
  }
  
  std::cout << " ------------ Ending tests ------------- " << std::endl;

  // ---------------------
  // Closing all the trees 
  // ---------------------
  treeCloudsA->close();
  treeCloudsB->close();
  treeHBookAout->close();
  treeHBookBout->close();

  std::cout << " END " << std::endl;

}
  
