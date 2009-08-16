{ 
  // ---------------------------------------------------------------------------
  // Last update: 13-Aug-2009.
  // 
  // This is the main ROOT (unnamed) macro.
  // It assumes that in the same directory there are 2 ROOT files:
  //   1)  ntuple_a.root
  //   2)  ntuple_b.root
  // each containing at least the tree "SimplifiedCalorimeter",
  // created by running the Geant4 application StatAccepTest. 
  // This macro loads the file  extraRootFunctions.C , where 
  // various functions are defined and then called.
  // It produces in output PostScript files for each observable
  // whose pvalue is below a threshold 
  // (see  pvalueThreshold ).
  //
  // ---------------------------------------------------------------------------

  // --------------------
  // Load what is needed.
  // --------------------
  // In order to use "G4double" instead of simply "double" we need
  // to tell Cint where to find the Geant4 include file G4Types.hh
  // and then we need to load the shared library of the dictionary.
  gInterpreter->AddIncludePath( "$G4INSTALL/source/global/management/include" );

  #include <vector>
  #include "G4Types.hh"

  gSystem->Load( "libG4TypesDict.so" );

  // Load the Root functions used by this script.
  gROOT->ProcessLine( ".L extraRootFunctions.C" );

  // ------------------------
  // Consider the first tree.
  // ------------------------
  TFile *file_a = new TFile( "ntuple_a.root" );
  TTree *tree_a = (TTree*) file_a->Get( "SimplifiedCalorimeter" );

  // For each variable in the ntuple we need to create a vector
  // whose element "i" is the value of that variable at the
  // event "i".
  // In the case of the longitudinal and transverse profile,
  // we would need a vector of a vector, but ROOT has problems
  // with that. We overcome the problem by using an array of
  // vectors, whose dimension is known by looking for instance
  // at the first event.

  std::vector< Float_t > vec_edep_act_a;
  std::vector< Float_t > vec_edep_cal_a;

  Int_t nLayers_a = 0;
  Int_t nBinR_a = 0;

  tree_a->SetBranchAddress( "nLayers", &nLayers_a );
  tree_a->SetBranchAddress( "nBinR", &nBinR_a );

  tree_a->GetEntry( 0 );

  std::vector< G4double > *pvecL_a[ nLayers_a ];
  for ( Int_t j = 0 ; j < nLayers_a ; j++ ) {
    pvecL_a[ j ] = new std::vector< G4double >;
  }

  std::vector< G4double > *pvecR_a[ nBinR_a ];
  for ( Int_t j = 0 ; j < nBinR_a ; j++ ) {
    pvecR_a[ j ] = new std::vector< G4double >;
  }

  getVectorsFromTree( tree_a, nLayers_a, nBinR_a,                           // input
		      vec_edep_act_a, vec_edep_cal_a, pvecL_a, pvecR_a );   // output

  // -------------------------
  // Consider the second tree.
  // -------------------------
  // Do exactly the same as for the first tree.
  TFile *file_b = new TFile( "ntuple_b.root" );
  TTree *tree_b = (TTree*) file_b->Get( "SimplifiedCalorimeter" );
  std::vector< Float_t > vec_edep_act_b;
  std::vector< Float_t > vec_edep_cal_b;
  Int_t nLayers_b = 0;
  Int_t nBinR_b = 0;
  tree_b->SetBranchAddress( "nLayers", &nLayers_b );
  tree_b->SetBranchAddress( "nBinR", &nBinR_b );
  tree_b->GetEntry( 0 );
  std::vector< G4double > *pvecL_b[ nLayers_b ];
  for ( Int_t j = 0 ; j < nLayers_b ; j++ ) {
    pvecL_b[ j ] = new std::vector< G4double >;
  }
  std::vector< G4double > *pvecR_b[ nBinR_b ];
  for ( Int_t j = 0 ; j < nBinR_b ; j++ ) {
    pvecR_b[ j ] = new std::vector< G4double >;
  }
  getVectorsFromTree( tree_b, nLayers_b, nBinR_b,                           // input
		      vec_edep_act_b, vec_edep_cal_b, pvecL_b, pvecR_b );   // output

  // -------------------------
  // Do the statistical tests.
  // -------------------------
  Double_t pvalue_edep_act = doStatisticalTest( vec_edep_act_a, vec_edep_act_b );

  cout << " pvalue_edep_act = " << pvalue_edep_act << endl;

  Double_t pvalue_edep_cal = doStatisticalTest( vec_edep_cal_a, vec_edep_cal_b );

  cout << " pvalue_edep_cal = " << pvalue_edep_cal << endl;

  std::vector< Double_t > pvalue_vecL( nLayers_a, -1.0 );
  if ( nLayers_a == nLayers_b ) {
    for ( Int_t j = 0 ; j < nLayers_a ; j++ ) {
      pvalue_vecL[ j ] = doStatisticalTest( *pvecL_a[ j ], *pvecL_b[ j ] );

      cout << " pvalue_vecL[" << j << "] = " << pvalue_vecL[ j ] << endl;
    }
  } else {
    cout << " ***ERROR: nLayers_a = " << nLayers_a 
	 << " NOT EQUAL TO nLayers_b = " << nLayers_b << endl;
  }
  
  std::vector< Double_t > pvalue_vecR( nBinR_a, -1.0 );
  if ( nBinR_a == nBinR_b ) {
    for ( Int_t j = 0 ; j < nBinR_a ; j++ ) {
      pvalue_vecR[ j ] = doStatisticalTest( *pvecR_a[ j ], *pvecR_b[ j ] );

      cout << " pvalue_vecR[" << j << "] = " << pvalue_vecR[ j ] << endl;
    }
  } else {
    cout << " ***ERROR: nBinR_a = " << nBinR_a 
	 << " NOT EQUAL TO nBinR_b = " << nBinR_b << endl;
  }

  // ---------------------------------------------------------
  // Produce a .ps file for each test below a given threshold.
  // ---------------------------------------------------------
  const Double_t pvalueThreshold = 0.01;   //***LOOKHERE***

  if ( pvalue_edep_act >= 0.0  &&  pvalue_edep_act < pvalueThreshold ) {
    doPlot( "edep_act", vec_edep_act_a, vec_edep_act_b );
  }

  if ( pvalue_edep_cal >= 0.0  &&  pvalue_edep_cal < pvalueThreshold ) {
    doPlot( "edep_cal", vec_edep_cal_a, vec_edep_cal_b );
  }

  for ( Int_t j = 0 ; j < nLayers_a ; j++ ) {
    if ( pvalue_vecL[ j ] >= 0.0  &&  pvalue_vecL[ j ] < pvalueThreshold ) {
      char buf[ 3 ];
      sprintf( buf, "L%i", j );
      doPlot( buf, *pvecL_a[ j ], *pvecL_b[ j ] );
    }
  }
  
  for ( Int_t j = 0 ; j < nBinR_a ; j++ ) {
    if ( pvalue_vecR[ j ] >= 0.0  &&  pvalue_vecR[ j ] < pvalueThreshold ) {
      char buf[ 3 ];
      sprintf( buf, "R%i", j );
      doPlot( buf, *pvecR_a[ j ], *pvecR_b[ j ] );
    }
  }

  // ---------------------
  // Delete, clean, close.
  // ---------------------
  delete [] pvecL_a;
  delete [] pvecR_a;
  delete [] pvecL_b;
  delete [] pvecR_b;

} // End unnamed ROOT macro.

