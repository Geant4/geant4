//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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
#include <vector>
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
 
public:
  //private constructor in order to create a singleton
 

  FCALAnalysisManager();
 
  G4String outputFileName;

  //  G4double OutOfWorld, Secondary, EmEdep, HadEdep; 

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

  AIDA::IHistogram1D* getfhisto_1() {return histo_1;} 
  AIDA::IHistogram1D* getfhisto_2() {return histo_2;} 
  AIDA::IHistogram1D* getfhisto_3() {return histo_3;} 
  AIDA::IHistogram1D* getfhisto_4() {return histo_4;} 

  AIDA::IHistogram1D*   histo_1;
  AIDA::IHistogram1D*   histo_2;
  AIDA::IHistogram1D*   histo_3;
  AIDA::IHistogram1D*   histo_4;

#endif
};

#endif



