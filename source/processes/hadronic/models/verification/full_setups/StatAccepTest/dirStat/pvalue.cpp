//========================================================================
//                         Program  pvalue.cpp
//                         -------------------
// This program reads in input two HBOOK files: 
//          1)  ntuple_a.hbook 
//          2)  ntuple_b.hbook
// which contain each a ntuple, containing physics information about
// the showering of a certain incoming beam particle into a sampling 
// calorimeter, and two histograms (that summarize, in a single plot,
// the longitudinal and transverse shower profile). 
// One of these two HBOOK files can be obtained with a certain version
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
//              from the first file (ntuple_a.hbook), and those obtained
//              from the second file (ntuple_b.hbook), respectively.
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
#include "CramerVonMisesBinnedComparisonAlgorithm.h"  
#include "AndersonDarlingBinnedComparisonAlgorithm.h" 

#include "KolmogorovSmirnovComparisonAlgorithm.h"         
#include "CramerVonMisesUnbinnedComparisonAlgorithm.h" 
#include "AndersonDarlingUnbinnedComparisonAlgorithm.h"

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
 
  // Define some useful variables.
  bool readOnly, createNew;
  int i_nLayers, i_nBinR, numberOfReplicas, numberOfRadiusBins; 
  int iEdepAct, iEdepCal, iL, iR; 

  // --- Start now working on the first ntuple ---

  // Read the first tree from the input file.
  readOnly  = true;  
  createNew = false; 
  std::auto_ptr<AIDA::ITree> treeHBookA( tf->create( "ntuple_a.hbook", "hbook", 
 						     readOnly, createNew ) );
  if ( ! treeHBookA.get() ) {
    std::cout << " ERROR: unable to open file  ntuple_a.hbook " << std::endl;
    return -1;
  } else {
    // std:: cout << " OK : opened file  ntuple_a.hbook " << std::endl   //***DEBUG***
    //           << "\t The content is: " << treeHBookA->ls() << std::endl;
  }
 
  // std::cout << " treeHBookA " << (int) treeHBookA.get() << std::endl; //***DEBUG***
  AIDA::ITuple* p_tpA = dynamic_cast<AIDA::ITuple*>( ( treeHBookA->find( "/1" ) ) );
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
 
  // Read also the two histograms: longitudinal shower profile
  // and transverse shower profile.
  AIDA::IHistogram1D& histL_A = 
    dynamic_cast<AIDA::IHistogram1D&>( * ( treeHBookA->find( "/50" ) ) );
  std::cout << " HistoL A title : " << histL_A.title() << std::endl;   //***DEBUG***
  AIDA::IHistogram1D& histR_A = 
    dynamic_cast<AIDA::IHistogram1D&>( * ( treeHBookA->find( "/60" ) ) );
  std::cout << " HistoR A title : " << histR_A.title() << std::endl;   //***DEBUG***
 
  // ---------------------------------------------------------------
  // Get the number of elements of the vector L and R of the ntuple.
  // ---------------------------------------------------------------  
  i_nLayers = tpA.findColumn( "nLayers" ); 
  i_nBinR = tpA.findColumn( "nBinR" ); 
  tpA.start(); tpA.next();
  numberOfReplicas   = tpA.getInt( i_nLayers ); 
  numberOfRadiusBins = tpA.getInt( i_nBinR ); 
  std::cout << " numberOfReplicas   = " << numberOfReplicas << std::endl
	    << " numberOfRadiusBins = " << numberOfRadiusBins 
 	    << std::endl;                                       //***DEBUG***
 
  // -------------------------------------------------------------
  // Creating the clouds, and fill them by projecting the ntuples.
  // -------------------------------------------------------------
  readOnly  = false;  
  createNew = true; 
  std::auto_ptr<AIDA::ITree> treeCloudsA( tf->create("cloudsA.xml", "xml", 
 						     readOnly, createNew) );
  
  // Creating a histogram factory attached to the tree.
  std::auto_ptr<AIDA::IHistogramFactory> 
    hfA( af->createHistogramFactory( *treeCloudsA ) );
  
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
  
  iEdepAct = tpA.findColumn( "EDEP_ACT" ); 
  iEdepCal = tpA.findColumn( "EDEP_CAL" ); 
  iL = tpA.findColumn( "L" );
  iR = tpA.findColumn( "R" ); 
   
  // Fill the clouds with the first ntuple.
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
  
  // Flushing the clouds into the file.
  treeCloudsA->commit();
  
  // --- Do the same for the ntuple B ---
  readOnly  = true;  
  createNew = false; 
  std::auto_ptr<AIDA::ITree> treeHBookB( tf->create( "ntuple_b.hbook", "hbook", 
                                                     readOnly, createNew ) );
  if ( ! treeHBookB.get() ) {
    std::cout << " ERROR: unable to open file  ntuple_b.hbook " << std::endl;
    return -1;
  }
  AIDA::ITuple& tpB = dynamic_cast<AIDA::ITuple&>( * ( treeHBookB->find( "/1" ) ) );
  AIDA::IHistogram1D& histL_B = 
    dynamic_cast<AIDA::IHistogram1D&>( * ( treeHBookB->find( "/50" ) ) );
  std::cout << " HistoL B title : " << histL_B.title() << std::endl;   //***DEBUG***
  AIDA::IHistogram1D& histR_B = 
    dynamic_cast<AIDA::IHistogram1D&>( * ( treeHBookB->find( "/60" ) ) );
  std::cout << " HistoR B title : " << histR_B.title() << std::endl;   //***DEBUG***
  i_nLayers = tpB.findColumn( "nLayers" ); 
  i_nBinR = tpB.findColumn( "nBinR" ); 
  tpB.start(); tpB.next();
  numberOfReplicas   = tpB.getInt( i_nLayers ); 
  numberOfRadiusBins = tpB.getInt( i_nBinR ); 
  readOnly  = false;  
  createNew = true; 
  std::auto_ptr<AIDA::ITree> treeCloudsB( tf->create("cloudsB.xml", "xml", 
  						     readOnly, createNew) );
  std::auto_ptr<AIDA::IHistogramFactory> 
    hfB( af->createHistogramFactory( *treeCloudsB ) );
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
  iEdepAct = tpB.findColumn( "EDEP_ACT" ); 
  iEdepCal = tpB.findColumn( "EDEP_CAL" ); 
  iL = tpB.findColumn( "L" );
  iR = tpB.findColumn( "R" ); 
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
  treeCloudsB->commit();
 
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
  AIDA::ICloud1D * cA2bis = 
    hfA->createCopy( "Energy deposit in All layers, A bis", *cA2 ); 
  std::vector< AIDA::ICloud1D * > cALbis;
  std::vector< AIDA::ICloud1D * > cARbis;
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    std::string name = "LA" + toString( i ) + "bis";
    cALbis.push_back( hfA->createCopy( name, *cAL[i] ) ); 
    // std::cout << "\t Created cloud: name = " << name << std::endl;  //***DEBUG*** 
  }
  for ( int i = 0; i < numberOfRadiusBins; i++ ) {
    std::string name = "RA" + toString( i ) + "bis";
    cARbis.push_back( hfA->createCopy( name, *cAR[i] ) ); 
    // std::cout << "\t Created cloud: name = " << name << std::endl;  //***DEBUG*** 
  }

  AIDA::ICloud1D * cB1bis = 
    hfB->createCopy( "Energy deposit in Active layers, B bis", *cB1 );
  AIDA::ICloud1D * cB2bis = 
    hfB->createCopy( "Energy deposit in All layers, B bis", *cB2 ); 
  std::vector< AIDA::ICloud1D * > cBLbis;
  std::vector< AIDA::ICloud1D * > cBRbis;
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    std::string name = "LB" + toString( i ) + "bis";
    cBLbis.push_back( hfB->createCopy( name, *cBL[i] ) ); 
    // std::cout << "\t Created cloud: name = " << name << std::endl;  //***DEBUG***
  }
  for ( int i = 0; i < numberOfRadiusBins; i++ ) {
    std::string name = "RB" + toString( i ) + "bis";
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
  // Creating a histogram factory, whose histograms will be handled by the tree
  std::auto_ptr<AIDA::IHistogramFactory> 
    hfacA( af->createHistogramFactory( *treeHBookAout ) );

  hfacA->createCopy( "11", cA1bis->histogram() );
  hfacA->createCopy( "12", cA2bis->histogram() );
  hfacA->createCopy( "50", histL_A );
  hfacA->createCopy( "60", histR_A );
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    int idNum = i + 100;
    std::string idString = toString( idNum );
    hfacA->createCopy( idString, cALbis[i]->histogram() );
  }
  for ( int i = 0; i < numberOfRadiusBins; i++ ) {
    int idNum = i + 200;
    std::string idString = toString( idNum );
    hfacA->createCopy( idString, cARbis[i]->histogram() );
  }
  
  // Flushing the histograms into the file
  treeHBookAout->commit();

  // --- Do the same for the ntuple B --- 
  std::auto_ptr<AIDA::ITree> treeHBookBout( tf->create("histo_b.hbook", "hbook", 
						       readOnly, createNew ) );   
  std::auto_ptr<AIDA::IHistogramFactory> 
    hfacB( af->createHistogramFactory( *treeHBookBout ) );
  hfacB->createCopy( "11", cB1bis->histogram() );
  hfacB->createCopy( "12", cB2bis->histogram() );
  hfacB->createCopy( "50", histL_B );
  hfacB->createCopy( "60", histR_B );
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    int idNum = i + 100;
    std::string idString = toString( idNum );
    hfacB->createCopy( idString, cBLbis[i]->histogram() );
  }
  for ( int i = 0; i < numberOfRadiusBins; i++ ) {
    int idNum = i + 200;
    std::string idString = toString( idNum );
    hfacB->createCopy( idString, cBRbis[i]->histogram() );
  }
  treeHBookBout->commit();

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

  // ------------------------
  // Do the statistical tests
  // ------------------------

  StatisticsComparator< Chi2ComparisonAlgorithm > comparatorBinC2;                    
  // StatisticsComparator< CramerVonMisesBinnedComparisonAlgorithm > comparatorBinCVM;
  // StatisticsComparator< AndersonDarlingBinnedComparisonAlgorithm > comparatorBinAD;

  StatisticsComparator< KolmogorovSmirnovComparisonAlgorithm > comparatorUnbinKS;
  StatisticsComparator< CramerVonMisesUnbinnedComparisonAlgorithm > comparatorUnbinCVM;
  StatisticsComparator< AndersonDarlingUnbinnedComparisonAlgorithm > comparatorUnbinAD;
  
  std::cout << " ------------ Starting tests ------------- " << std::endl;
  
  std::cout << "Observable 1";
  ComparisonResult resultBinC2_1 = comparatorBinC2.compare( cA1bis->histogram(), 
							    cB1bis->histogram() );
  ComparisonResult resultUnbinKS_1  = comparatorUnbinKS.compare( *cA1, *cB1 ); 
  ComparisonResult resultUnbinCVM_1 = comparatorUnbinCVM.compare( *cA1, *cB1 ); 
  ComparisonResult resultUnbinAD_1  = comparatorUnbinAD.compare( *cA1, *cB1 ); 
  if ( resultBinC2_1.quality()    < pvalueThreshold  || 
       resultUnbinKS_1.quality()  < pvalueThreshold  ||
       resultUnbinCVM_1.quality() < pvalueThreshold  ||
       resultUnbinAD_1.quality()  < pvalueThreshold ) {
    std::cout << "\t ***WARNING***";
  }
  std::cout << std::endl;
  std::cout << "  Chi2" << "  d=" << resultBinC2_1.distance()
	    << "  ndf=" << resultBinC2_1.ndf() 
	    << "  pvalue=" << resultBinC2_1.quality() << std::endl;
  std::cout << "  KS" << "  d=" << resultUnbinKS_1.distance()
	    << "  ndf=" << resultUnbinKS_1.ndf() 
	    << "  pvalue=" << resultUnbinKS_1.quality() << std::endl;
  std::cout << "  CVM" << "  d=" << resultUnbinCVM_1.distance()
	    << "  ndf=" << resultUnbinCVM_1.ndf() 
	    << "  pvalue=" << resultUnbinCVM_1.quality() << std::endl;
  std::cout << "  AD" << "  d=" << resultUnbinAD_1.distance()
	    << "  ndf=" << resultUnbinAD_1.ndf() 
	    << "  pvalue=" << resultUnbinAD_1.quality() << std::endl;
      
  std::cout << "Observable 2";
  ComparisonResult resultBinC2_2 = comparatorBinC2.compare( cA2bis->histogram(), 
							    cB2bis->histogram() );
  ComparisonResult resultUnbinKS_2  = comparatorUnbinKS.compare( *cA2, *cB2 ); 
  ComparisonResult resultUnbinCVM_2 = comparatorUnbinCVM.compare( *cA2, *cB2 ); 
  ComparisonResult resultUnbinAD_2  = comparatorUnbinAD.compare( *cA2, *cB2 ); 
  if ( resultBinC2_2.quality()    < pvalueThreshold  || 
       resultUnbinKS_2.quality()  < pvalueThreshold  ||
       resultUnbinCVM_2.quality() < pvalueThreshold  ||
       resultUnbinAD_2.quality()  < pvalueThreshold ) {
    std::cout << "\t ***WARNING***";
  }
  std::cout << std::endl;
  std::cout << "  Chi2" << "  d=" << resultBinC2_2.distance()
	    << "  ndf=" << resultBinC2_2.ndf() 
	    << "  pvalue=" << resultBinC2_2.quality() << std::endl;
  std::cout << "  KS" << "  d=" << resultUnbinKS_2.distance()
            << "  ndf=" << resultUnbinKS_2.ndf() 
            << "  pvalue=" << resultUnbinKS_2.quality() << std::endl;
  std::cout << "  CVM" << "  d=" << resultUnbinCVM_2.distance()
            << "  ndf=" << resultUnbinCVM_2.ndf() 
            << "  pvalue=" << resultUnbinCVM_2.quality() << std::endl;
  std::cout << "  AD" << "  d=" << resultUnbinAD_2.distance()
            << "  ndf=" << resultUnbinAD_2.ndf() 
            << "  pvalue=" << resultUnbinAD_2.quality() << std::endl;

  // For the time being, we comment out the tests for the longitudinal 
  // and transverse shower shapes because the error bars are not correct, 
  // due to a "feature" of PI which does not allow us to set the proper 
  // errors.

  // std::cout << "Observable 3";  // Longitudinal shower shape: Only binned tests.
  // ComparisonResult resultBinC2_3  = comparatorBinC2.compare( histL_A, histL_B );
  // ComparisonResult resultBinCVM_3 = comparatorBinCVM.compare( histL_A, histL_B );
  // ComparisonResult resultBinAD_3  = comparatorBinAD.compare( histL_A, histL_B );
  // if ( resultBinC2_3.quality()  < pvalueThreshold  ||
  //      resultBinCVM_3.quality() < pvalueThreshold  ||
  //      resultBinAD_3.quality()  < pvalueThreshold ) {
  //   std::cout << "\t ***WARNING***";
  // }
  // std::cout << std::endl;
  // std::cout << "  Chi2" << "  d=" << resultBinC2_3.distance()
  // 	       << "  ndf=" << resultBinC2_3.ndf() 
  // 	       << "  pvalue=" << resultBinC2_3.quality() << std::endl;
  // std::cout << "  CVM binned" << "  d=" << resultBinCVM_3.distance()
  // 	       << "  ndf=" << resultBinCVM_3.ndf() 
  // 	       << "  pvalue=" << resultBinCVM_3.quality() << std::endl;
  // std::cout << "  AD binned" << "  d=" << resultBinAD_3.distance()
  // 	       << "  ndf=" << resultBinAD_3.ndf() 
  // 	       << "  pvalue=" << resultBinAD_3.quality() << std::endl;
  // 
  // std::cout << "Observable 4";  // Transverse shower shape: Only binned tests.
  // ComparisonResult resultBinC2_4  = comparatorBinC2.compare( histR_A, histR_B );
  // ComparisonResult resultBinCVM_4 = comparatorBinCVM.compare( histR_A, histR_B );
  // ComparisonResult resultBinAD_4  = comparatorBinAD.compare( histR_A, histR_B );
  // if ( resultBinC2_4.quality()  < pvalueThreshold  ||
  //      resultBinCVM_4.quality() < pvalueThreshold  ||
  //      resultBinAD_4.quality()  < pvalueThreshold ) {
  //   std::cout << "\t ***WARNING***";
  // }
  // std::cout << std::endl;
  // std::cout << "  Chi2" << "  d=" << resultBinC2_4.distance()
  // 	       << "  ndf=" << resultBinC2_4.ndf() 
  // 	       << "  pvalue=" << resultBinC2_4.quality() << std::endl;
  // std::cout << "  CVM binned" << "  d=" << resultBinCVM_4.distance()
  // 	       << "  ndf=" << resultBinCVM_4.ndf() 
  // 	       << "  pvalue=" << resultBinCVM_4.quality() << std::endl;
  // std::cout << "  AD binned" << "  d=" << resultBinAD_4.distance()
  // 	       << "  ndf=" << resultBinAD_4.ndf() 
  // 	       << "  pvalue=" << resultBinAD_4.quality() << std::endl;
  
  std::vector< ComparisonResult > resultBinC2_L;
  std::vector< ComparisonResult > resultUnbinKS_L;
  std::vector< ComparisonResult > resultUnbinCVM_L;
  std::vector< ComparisonResult > resultUnbinAD_L;
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    std::cout << "Observable L" << i;   
    ComparisonResult resultBinC2_Lvalue = 
      comparatorBinC2.compare( cALbis[i]->histogram(), cBLbis[i]->histogram() );
    resultBinC2_L.push_back( resultBinC2_Lvalue );

    ComparisonResult resultUnbinKS_Lvalue = 
      comparatorUnbinKS.compare( *cAL[i], *cBL[i] );
    resultUnbinKS_L.push_back( resultUnbinKS_Lvalue );
    ComparisonResult resultUnbinCVM_Lvalue = 
      comparatorUnbinCVM.compare( *cAL[i], *cBL[i] );
    resultUnbinCVM_L.push_back( resultUnbinCVM_Lvalue );
    ComparisonResult resultUnbinAD_Lvalue = 
      comparatorUnbinAD.compare( *cAL[i], *cBL[i] );
    resultUnbinAD_L.push_back( resultUnbinAD_Lvalue );
    if ( resultBinC2_Lvalue.quality()    < pvalueThreshold ||
	 resultUnbinKS_Lvalue.quality()  < pvalueThreshold ||
	 resultUnbinCVM_Lvalue.quality() < pvalueThreshold ||
	 resultUnbinAD_Lvalue.quality()  < pvalueThreshold ) {
      std::cout << "\t ***WARNING***";
    }
    std::cout << std::endl;
    std::cout << "  Chi2" << "  d=" << resultBinC2_Lvalue.distance()
  	      << "  ndf=" << resultBinC2_Lvalue.ndf() 
  	      << "  pvalue=" << resultBinC2_Lvalue.quality() << std::endl;
    std::cout << "  KS" << "  d=" << resultUnbinKS_Lvalue.distance()
  	      << "  ndf=" << resultUnbinKS_Lvalue.ndf() 
  	      << "  pvalue=" << resultUnbinKS_Lvalue.quality() << std::endl;
    std::cout << "  CVM" << "  d=" << resultUnbinCVM_Lvalue.distance()
  	      << "  ndf=" << resultUnbinCVM_Lvalue.ndf() 
  	      << "  pvalue=" << resultUnbinCVM_Lvalue.quality() << std::endl;
    std::cout << "  AD" << "  d=" << resultUnbinAD_Lvalue.distance()
  	      << "  ndf=" << resultUnbinAD_Lvalue.ndf() 
  	      << "  pvalue=" << resultUnbinAD_Lvalue.quality() << std::endl;
  }

  std::vector< ComparisonResult > resultBinC2_R;
  std::vector< ComparisonResult > resultUnbinKS_R;
  std::vector< ComparisonResult > resultUnbinCVM_R;
  std::vector< ComparisonResult > resultUnbinAD_R;
  for ( int i = 0; i < numberOfRadiusBins; i++ ) {
    std::cout << "Observable R" << i;   
    ComparisonResult resultBinC2_Rvalue = 
      comparatorBinC2.compare( cARbis[i]->histogram(), cBRbis[i]->histogram() );
    resultBinC2_R.push_back( resultBinC2_Rvalue );

    ComparisonResult resultUnbinKS_Rvalue = 
      comparatorUnbinKS.compare( *cAR[i], *cBR[i] );
    resultUnbinKS_R.push_back( resultUnbinKS_Rvalue );
    ComparisonResult resultUnbinCVM_Rvalue = 
      comparatorUnbinCVM.compare( *cAR[i], *cBR[i] );
    resultUnbinCVM_R.push_back( resultUnbinCVM_Rvalue );
    ComparisonResult resultUnbinAD_Rvalue = 
      comparatorUnbinAD.compare( *cAR[i], *cBR[i] );
    resultUnbinAD_R.push_back( resultUnbinAD_Rvalue );
    if ( resultBinC2_Rvalue.quality()    < pvalueThreshold ||
	 resultUnbinKS_Rvalue.quality()  < pvalueThreshold ||
	 resultUnbinCVM_Rvalue.quality() < pvalueThreshold ||
	 resultUnbinAD_Rvalue.quality()  < pvalueThreshold ) {
      std::cout << "\t ***WARNING***";
    }
    std::cout << std::endl;
    std::cout << "  Chi2" << "  d=" << resultBinC2_Rvalue.distance()
  	      << "  ndf=" << resultBinC2_Rvalue.ndf() 
  	      << "  pvalue=" << resultBinC2_Rvalue.quality() << std::endl;
    std::cout << "  KS" << "  d=" << resultUnbinKS_Rvalue.distance()
  	      << "  ndf=" << resultUnbinKS_Rvalue.ndf() 
  	      << "  pvalue=" << resultUnbinKS_Rvalue.quality() << std::endl;
    std::cout << "  CVM" << "  d=" << resultUnbinCVM_Rvalue.distance()
  	      << "  ndf=" << resultUnbinCVM_Rvalue.ndf() 
  	      << "  pvalue=" << resultUnbinCVM_Rvalue.quality() << std::endl;
    std::cout << "  AD" << "  d=" << resultUnbinAD_Rvalue.distance()
  	      << "  ndf=" << resultUnbinAD_Rvalue.ndf() 
  	      << "  pvalue=" << resultUnbinAD_Rvalue.quality() << std::endl;
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
  
