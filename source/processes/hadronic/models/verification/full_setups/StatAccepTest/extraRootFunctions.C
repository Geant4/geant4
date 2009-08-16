#include <algorithm>
#include <vector>
#include "G4Types.hh"
#include <TMath.h>
#include <TH1D.h>
#include <Rtypes.h>


// This function receives in input a tree, and two integers.
// It returns 2 vectors of floatings: 
//   - visible energy 
//   - total energy 
// for each event, ordered in ascending order;
// and 2 vectors of arrays of G4double:
//   - visible energy in each longitudinal layer
//   - visible energy in each transverse ring
// for each event, ordered in ascending order.
void getVectorsFromTree( const TTree *tree,                             // input
			 const Int_t nLayers,                           // input
			 const Int_t nBinR,                             // input
			 std::vector< Float_t > & vec_edep_act,         // output
			 std::vector< Float_t > & vec_edep_cal,         // output
			 std::vector< G4double > *pvecL[ nLayers ],     // output
			 std::vector< G4double > *pvecR[ nBinR ]        // output
		       ) {
  Int_t id = 0;
  Double_t e = 0.0;
  Float_t edep_act = 0.0;
  Float_t edep_cal = 0.0;

  std::vector< G4double > *plongitudinalProfile = new std::vector< G4double >;
  std::vector< G4double > *ptransverseProfile = new std::vector< G4double >;

  tree->SetBranchAddress( "ID", &id );
  tree->SetBranchAddress( "E", &e );
  tree->SetBranchAddress( "EDEP_ACT", &edep_act );
  tree->SetBranchAddress( "EDEP_CAL", &edep_cal );
  tree->SetBranchAddress( "L", &plongitudinalProfile );
  tree->SetBranchAddress( "R", &ptransverseProfile );
  
  Int_t numEntries = ( Int_t ) tree->GetEntries();

  for ( Int_t i = 0 ; i < numEntries ; i++ ) {  // Event loop
  
    tree->GetEntry( i );

    //cout << "Event i=" << i << endl                 //***DEBUG***
    //     << "\t id=" << id << endl
    //	 << "\t e=" << e << endl
    //     << "\t edep_act=" << edep_act << endl
    //     << "\t edep_cal=" << edep_cal << endl;
    //cout << "\t L : " << endl;
    //for ( Int_t j = 0 ; j < nLayers ; j++ ) { 
    //  cout << "\t \t L[" << j << "]=" << ( *plongitudinalProfile )[ j ] << endl;
    //}
    //cout << "\t R : " << endl;
    //for ( Int_t j = 0 ; j < nBinR ; j++ ) { 
    //  cout << "\t \t R[" << j << "]=" << ( *ptransverseProfile )[ j ] << endl;
    //}
    //cout << endl;

    vec_edep_act.push_back( edep_act );
    vec_edep_cal.push_back( edep_cal );

    for ( Int_t j = 0 ; j < nLayers ; j++ ) { 
      G4double val = ( *plongitudinalProfile )[ j ];
      pvecL[ j ]->push_back( val );
    }
    
    for ( Int_t j = 0 ; j < nBinR ; j++ ) { 
      G4double val = ( *ptransverseProfile )[ j ];
      pvecR[ j ]->push_back( val );
    }

  } // End event loop

  //***DEBUG*** : Check if the vectors are ok
  //cout << endl << " === CHECK VECTORS === " << endl << endl;
  //Int_t count = 0;
  //cout << endl << "  === vec_edep_act === " << endl;
  //for ( std::vector< Float_t >::const_iterator cit_act = vec_edep_act.begin() ;
  //	cit_act != vec_edep_act.end() ; ++cit_act ) {
  //  cout << "\t" << count++ << "\t" << *cit_act << endl;
  //}
  //count = 0;
  //cout << endl << "  === vec_edep_cal === " << endl;
  //for ( std::vector< Float_t >::const_iterator cit_cal = vec_edep_cal.begin() ;
  //	cit_cal != vec_edep_cal.end() ; ++cit_cal ) {
  //  cout << "\t" << count++ << "\t" << *cit_cal << endl;
  //}
  //for ( Int_t j = 0 ; j < nLayers ; j++ ) { 
  //  std::vector< G4double > v( ( *pvecL[ j ] ) );
  //  cout << endl << "  === vecL === " << j << endl;
  //  count = 0;
  //  for ( std::vector< G4double >::const_iterator cit = v.begin() ;
  //	  cit != v.end() ; ++cit ) {
  //    cout << "\t" << count++ << "\t" << *cit << endl;
  //  }
  //}
  //for ( Int_t j = 0 ; j < nBinR ; j++ ) { 
  //  std::vector< G4double > v( ( *pvecR[ j ] ) );
  //  cout << endl << "  === vecR === " << j << endl;
  //  count = 0;
  //  for ( std::vector< G4double >::const_iterator cit = v.begin() ;
  //	  cit != v.end() ; ++cit ) {
  //    cout << "\t" << count++ << "\t" << *cit << endl;
  //  }
  //}

  // We need to order the vectors (with event entries) in ascending order
  // to be able to use the unbinned Kolmogorov statistical test.
  
  std::sort( vec_edep_act.begin(), vec_edep_act.end() );
  std::sort( vec_edep_cal.begin(), vec_edep_cal.end() );
  for ( Int_t j = 0 ; j < nLayers ; j++ ) { 
    std::sort( pvecL[ j ]->begin(), pvecL[ j ]->end() );
  }
  for ( Int_t j = 0 ; j < nBinR ; j++ ) { 
    std::sort( pvecR[ j ]->begin(), pvecR[ j ]->end() );
  }

  //***DEBUG*** : Check if the vectors are ok after ordering
  //cout << endl << " === CHECK VECTORS AFTER ORDERING === " << endl << endl;
  //count = 0;
  //cout << endl << "  === vec_edep_act === " << endl;
  //for ( std::vector< Float_t >::const_iterator cit_act = vec_edep_act.begin() ;
  //	cit_act != vec_edep_act.end() ; ++cit_act ) {
  //  cout << "\t" << count++ << "\t" << *cit_act << endl;
  //}
  //count = 0;
  //cout << endl << "  === vec_edep_cal === " << endl;
  //for ( std::vector< Float_t >::const_iterator cit_cal = vec_edep_cal.begin() ;
  //	cit_cal != vec_edep_cal.end() ; ++cit_cal ) {
  //  cout << "\t" << count++ << "\t" << *cit_cal << endl;
  //}
  //for ( Int_t j = 0 ; j < nLayers ; j++ ) { 
  //  std::vector< G4double > v( ( *pvecL[ j ] ) );
  //  cout << endl << "  === vecL === " << j << endl;
  //  count = 0;
  //  for ( std::vector< G4double >::const_iterator cit = v.begin() ;
  //	  cit != v.end() ; ++cit ) {
  //    cout << "\t" << count++ << "\t" << *cit << endl;
  //  }
  //}
  //for ( Int_t j = 0 ; j < nBinR ; j++ ) { 
  //  std::vector< G4double > v( ( *pvecR[ j ] ) );
  //  cout << endl << "  === vecR === " << j << endl;
  //  count = 0;
  //  for ( std::vector< G4double >::const_iterator cit = v.begin() ;
  //	  cit != v.end() ; ++cit ) {
  //    cout << "\t" << count++ << "\t" << *cit << endl;
  //  }
  //}

  // Finally, release memory.
  delete plongitudinalProfile;
  delete ptransverseProfile;
}


Double_t doStatisticalTest( std::vector< Float_t > vec_a ,
			    std::vector< Float_t > vec_b ) {
  // Convert the vectors of Float_t into vectors of Double_t. 
  std::vector< Double_t > dvec_a;
  for ( std::vector< Float_t >::const_iterator cit = vec_a.begin();
	cit != vec_a.end(); ++cit ) {
    dvec_a.push_back( Double_t ( *cit ) );
  }
  std::vector< Double_t > dvec_b;
  for ( std::vector< Float_t >::const_iterator cit = vec_b.begin();
	cit != vec_b.end(); ++cit ) {
    dvec_b.push_back( Double_t ( *cit ) );
  }
  return doStatisticalTest( dvec_a , dvec_b );
}


Double_t doStatisticalTest( std::vector< G4double > vec_a ,
			    std::vector< G4double > vec_b ) {
  // Convert the vectors of G4double into vectors of Double_t. 
  std::vector< Double_t > dvec_a;
  for ( std::vector< G4double >::const_iterator cit = vec_a.begin();
	cit != vec_a.end(); ++cit ) {
    dvec_a.push_back( Double_t ( *cit ) );
  }
  std::vector< Double_t > dvec_b;
  for ( std::vector< G4double >::const_iterator cit = vec_b.begin();
	cit != vec_b.end(); ++cit ) {
    dvec_b.push_back( Double_t ( *cit ) );
  }
  return doStatisticalTest( dvec_a , dvec_b );
}


Double_t doStatisticalTest( std::vector< Double_t > vec_a ,
			    std::vector< Double_t > vec_b ) {

  // This function returns the lowest pvalue between a set
  // of statistical tests on the input vectors vec_a and vec_b.
  // We assume that the vectors vec_a and vec_b are already
  // ordered in ascending order.
  // For the time being, we only use the Kolmogorov test
  // (perhaps in the future we could add the Cramer-von Mises
  //  and the Anderson-Darling tests).

  Double_t pvalue = -1.0;

  // To call TMath::Kolmogorov the vectors should be converted to arrays.
  Double_t* pa = new Double_t[ vec_a.size() ]; 
  Int_t i = 0;
  for ( std::vector< Double_t >::const_iterator cit = vec_a.begin();
	cit != vec_a.end(); ++cit ) {
    pa[ i++ ] = *cit;
  }
  Double_t* pb = new Double_t[ vec_b.size() ]; 
  i = 0;
  for ( std::vector< Double_t >::const_iterator cit = vec_b.begin();
	cit != vec_b.end(); ++cit ) {
    pb[ i++ ] = *cit;
  }

  pvalue = TMath::KolmogorovTest( vec_a.size(), pa, vec_b.size(), pb, "" );

  // Finally, release the memory.
  delete [] pa;
  delete [] pb;

  return pvalue;
}


void doPlot( char *str, std::vector< Float_t > vec_a , std::vector< Float_t > vec_b ) {
  // Convert the vectors of Float_t into vectors of Double_t. 
  std::vector< Double_t > dvec_a;
  for ( std::vector< Float_t >::const_iterator cit = vec_a.begin();
	cit != vec_a.end(); ++cit ) {
    dvec_a.push_back( Double_t ( *cit ) );
  }
  std::vector< Double_t > dvec_b;
  for ( std::vector< Float_t >::const_iterator cit = vec_b.begin();
	cit != vec_b.end(); ++cit ) {
    dvec_b.push_back( Double_t ( *cit ) );
  }
  return doPlot( str, dvec_a , dvec_b );
}


void doPlot( char *str, std::vector< G4double > vec_a , std::vector< G4double > vec_b ) {
  // Convert the vectors of G4double into vectors of Double_t. 
  std::vector< Double_t > dvec_a;
  for ( std::vector< G4double >::const_iterator cit = vec_a.begin();
	cit != vec_a.end(); ++cit ) {
    dvec_a.push_back( Double_t ( *cit ) );
  }
  std::vector< Double_t > dvec_b;
  for ( std::vector< G4double >::const_iterator cit = vec_b.begin();
	cit != vec_b.end(); ++cit ) {
    dvec_b.push_back( Double_t ( *cit ) );
  }
  return doPlot( str, dvec_a , dvec_b );
}


void doPlot( char *str, std::vector< Double_t > vec_a , std::vector< Double_t > vec_b ) {

  // This function receives in input a string, and two vectors of
  // doubles (assumed to be ordered in ascending orders.
  // It produces a  .ps  file where two histograms are compared:
  //   RED  : corresponding to the first vector,  vec_a;
  //   BLUE : corresponding to the second vector, vec_b.
  // The histograms are obtained by dividing in equal-size bins
  // the interval between the lowest and largest value of the two
  // vectors (actually, an interval a bit wider than that is used,
  // to make sure that all events fit inside the plot).
  // The input string is used as title of the histogram.

  const Int_t Nbin = 100;  //***LOOKHERE*** : Number of bins for plotting.

  Double_t infX = vec_a[ 0 ];
  if ( vec_b[ 0 ] < infX ) infX = vec_b[ 0 ];
  infX = infX * 0.999 - 0.001;  // A bit smaller than the smallest value.

  Double_t supX = vec_a[ vec_a.size() - 1 ]; 
  if ( vec_b[ vec_b.size() - 1 ] > supX ) supX = vec_b[ vec_b.size() - 1 ];
  supX = supX * 1.001 + 0.001;  // A bit larger than the largest value.

  Double_t binSize = ( supX - infX ) / ( 1.0*Nbin );

  //cout << "Observable: " << str << endl                              //***DEBUG***
  //     << "\t #entries: " << vec_a.size() << "\t" << vec_b.size() << endl
  //     << "\t vec[0]:   " << vec_a[ 0 ]   << "\t" << vec_b[ 0 ] << endl
  //     << "\t vec[N-1]:   " << vec_a[ vec_a.size() - 1 ]   
  //     << "\t" << vec_b[ vec_b.size() - 1 ] << endl
  //     << "\t range: ( " << infX << " , " << supX << " ) "
  //     << "\t #bins: " << Nbin << "  size: " << binSize << endl;

  // Create the two histograms.
  char id[ 100 ];
  sprintf( id, "%s_%s", str, "a");
  //cout << " id=" << id << endl;
  TH1D *histo_a = new TH1D( id , str , Nbin , infX , supX );
  char title[ 100 ];
  sprintf( title, "%s %s", str, "(MeV)" );
  //cout << title << endl;
  histo_a->GetXaxis()->SetTitle( title );
  histo_a->GetYaxis()->SetTitle( "Events" );
 
  sprintf( id, "%s_%s", str, "b");  
  TH1D *histo_b = new TH1D( id , str , Nbin , infX , supX );
  histo_b->GetXaxis()->SetTitle( title );
  histo_b->GetYaxis()->SetTitle( "Events" );
  
  // Fill the two histograms.
  for ( std::vector< Double_t >::const_iterator cit = vec_a.begin();
  	cit != vec_a.end(); ++cit ) {
    histo_a->Fill( *cit );
  }
  for ( std::vector< Double_t >::const_iterator cit = vec_b.begin();
        cit != vec_b.end(); ++cit ) {
    histo_b->Fill( *cit );
  }
  
  // Plot and save it as a .ps file.
  TCanvas* c = new TCanvas( str, str );
  c->Divide( 1, 2 );
  c->cd( 1 );
  c->SetLogy( 0 );  // Linear scale

  gStyle->SetHistLineWidth( 3 );
  gStyle->SetLabelSize( 0.05 );
  gStyle->SetLabelSize( 0.05, "Y" );
  gStyle->SetOptStat( "e" );   // Only #entries
  histo_b->SetStats( kFALSE );
  gStyle->SetHistLineColor( kRed ); histo_a->UseCurrentStyle(); histo_a->Draw();
  gStyle->SetHistLineColor( kBlue ); histo_b->UseCurrentStyle(); histo_b->Draw( "same" );

  c->cd( 2 );
  c->SetLogy();     // Logarithmic scale
  gStyle->SetHistLineColor( kRed ); histo_a->UseCurrentStyle(); histo_a->Draw();
  gStyle->SetHistLineColor( kBlue ); histo_b->UseCurrentStyle(); histo_b->Draw( "same" );

  char plotFileName[ 100 ];
  sprintf( plotFileName, "%s_%s.%s", "plot", str, "ps" );
  //cout << plotFileName << endl;
  c->Print( plotFileName );
 
  // Finally, release the memory.
  delete histo_a;
  delete histo_b;
  delete c;

}
