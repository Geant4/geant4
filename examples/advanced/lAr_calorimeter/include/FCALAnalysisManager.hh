// $Id: FCALAnalysisManager.hh
// GEANT4 tag $Name: xray_fluo-V04-01-03 
//
// Author: Patricia Mendez (patricia.mendez@cern.ch)
//
// History:
// -----------
//  12 Feb 2003 Patricia Mendez created based on XrayFluoAnalysisManager.
// -------------------------------------------------------------------

#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH


#include "globals.hh"
#include "g4std/vector"
#include "G4ThreeVector.hh"
#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#endif
#include "FCALAnalysisMessenger.hh"

class G4Step;
#ifdef G4ANALYSIS_USE
namespace AIDA {
class IAnalysisFactory;
class IHistogramFactory;
class ITree;
class ITupleFactory;
class ITuple;
};
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FCALAnalysisManager
{
public:
 
  virtual ~FCALAnalysisManager();
  
#ifdef G4ANALYSIS_USE

  void book();
  
  void finish();
  
  
 //fill histograms with data from FCALTBEventAction
  void analyseEnergyDep(G4double eDep);

  //method to call to create an instance of this class
  static FCALAnalysisManager* getInstance();

  //method intended to chenge the name of the hbook output file
  void SetOutputFileName(G4String);
#endif
 
private:
  //private constructor in order to create a singleton
 

  FCALAnalysisManager();
 
  G4String outputFileName;

  G4double OutOfWorld, Secondary, EmEdep, HadEdep; 

  static FCALAnalysisManager* instance;

#ifdef G4ANALYSIS_USE  
  //pointer to the analysis messenger
  FCALAnalysisMessenger* analisysMessenger;

  AIDA::IAnalysisFactory* analysisFactory;
  AIDA::ITree* tree;

  AIDA::IHistogramFactory *histogramFactory;
  AIDA::  ITupleFactory* tupleFactory;

  AIDA::ITuple* ntuple_1;
  AIDA::ITuple* ntuple_2;
  AIDA::ITuple* ntuple_3;

  AIDA::IHistogram1D*   histo_1;
  AIDA::IHistogram1D*   histo_2;
  AIDA::IHistogram1D*   histo_3;
  AIDA::IHistogram1D*   histo_4;

#endif
};

#endif



