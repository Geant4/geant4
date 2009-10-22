//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DMXAnalysisManager.hh
// GEANT4 tag $Name:
//
// Author: Alex Howard (alexander.howard@cern.ch)
//
// History:
// -----------
//  16 Jan 2002  Alex Howard   Created
//  22 Oct 2009  Luciano Pandola Added TreeFactory, removed ntuple4
//
// -------------------------------------------------------------------



#ifdef G4ANALYSIS_USE
#ifndef DMXAnalysisManager_h
#define DMXAnalysisManager_h 1

#include "globals.hh"

// Histogramming from AIDA 


#include "AIDA/IAnalysisFactory.h"

#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"

#include "AIDA/IHistogramFactory.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
#include "AIDA/IHistogram3D.h"

// #include "AIDA/IPlotterFactory.h"
// #include "AIDA/IPlotterRegion.h"
// #include "AIDA/IPlotter.h"

#include "AIDA/ITupleFactory.h"
#include "AIDA/ITuple.h"

#include "AIDA/IManagedObject.h"

// # include "AIDA/IFitter.h"
// # include "AIDA/IFitResult.h"
// # include "AIDA/IFitData.h"
// # include "AIDA/IRangeSet.h"
// # include "AIDA/IFitParameterSettings.h"
// # include "AIDA/IFunctionFactory.h"
// # include "AIDA/IFunction.h"
// # include "AIDA/IFitFactory.h"

namespace AIDA {
  class IAnalysisFactory;
  class ITree;
  class ITreeFactory;
  class IHistogramFactory;
  class ITupleFactory;
  class ITuple;
  class IHistogram1D;
  class IHistogram2D;
  //class IPlotter;
  //  class IFitter;
  //class IFitResult;
  //class IFitData;
  //class IRangeSet;
  //class IFitParameterSettings;
  //class IFunctionFactory;
  //class IFunction;
  //class IFitFactory;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class DMXAnalysisManager
{
public:
 
  virtual ~DMXAnalysisManager();
  
  //  void book();
  void book(G4String, G4bool);
  
  void Init();
  void Finish();
  
  //fill histograms with SHC (Scint Hits) data from DMXEventAction
  void analyseScintHits(G4int event_id, G4double energy_pri, G4double totEnergy, G4int S_hits, G4double firstLXeHitTime, G4int P_hits, G4double aveTimePmtHits, G4String firstparticleName, G4double firstParticleE, G4bool gamma_ev, G4bool neutron_ev, G4bool positron_ev, G4bool electron_ev, G4bool other_ev, long seed1, long seed2);

  //fill histograms with PHC (Pmt Hits) data from DMXEventAction
  void analysePMTHits(G4int, G4int, G4double, G4double, G4double);
  
  //fill histograms with data from DMXParticleSource / secondary history
  void analyseParticleSource(G4double, G4String);

  //fill histograms with data from DMXPrimaryGenerator
  void analysePrimaryGenerator(G4double);
  
  //fill histograms with data from DMXPrimaryGenerator
  void HistTime(G4double);
  
  //fit exponential decay time of scintillation time:
  void PulseTimeFit();

  //Method to interactively display histograms
  void PlotHistos(G4bool);
  
  //Method to interactively display histograms
  void PlotHistosInit();
  //Method to interactively display histograms
  void PlotHistosInter(G4int flag);
  
  //method to call to create an instance of this class
  static DMXAnalysisManager* getInstance();
 

private:
  
  //private constructor in order to create a singleton
  DMXAnalysisManager();
 
  static DMXAnalysisManager* instance;
  
  // Quantities for the ntuple


  G4int event_id; 
  G4double energy_pri; 
  G4double totEnergy; 
  G4int S_hits; 
  G4double firstLXeHitTime; 
  G4int P_hits; 
  G4double aveTimePmtHits; 
  G4String firstparticleName; 
  G4double firstParticleE; 
  G4bool gamma_ev; 
  G4bool neutron_ev; 
  G4bool positron_ev;
  G4bool electron_ev;
  G4bool other_ev; 
  long seed1; long seed2;

  G4int i; G4double x; G4double y; G4double z;

  G4double energy; G4double name;

  G4double time;

  G4bool interactive;

  AIDA::IAnalysisFactory  *af;
  AIDA::ITreeFactory *tf;
  AIDA::ITree             *tree;
  AIDA::IHistogramFactory *hf;
  AIDA::ITupleFactory     *tpf;
  AIDA::ITuple            *ntuple1;
  AIDA::ITuple            *ntuple2;
  AIDA::ITuple            *ntuple3;


  AIDA::IHistogram1D* hEsourcep;
  AIDA::IHistogram1D* hEdepp;
  AIDA::IHistogram1D* hEdepRecoil;
  AIDA::IHistogram1D* hNumPhLow;
  AIDA::IHistogram1D* hNumPhHigh;
  AIDA::IHistogram1D* hAvPhArrival;
  AIDA::IHistogram1D* hPhArrival;
  AIDA::IHistogram2D* hPMTHits;
  AIDA::IHistogram2D* h1stPMTHit;
  AIDA::IHistogram1D* hGammaEdep;
  AIDA::IHistogram1D* hNeutronEdep;
  AIDA::IHistogram1D* hElectronEdep;
  AIDA::IHistogram1D* hPositronEdep;
  AIDA::IHistogram1D* hOtherEdep;

  //AIDA::IFunctionFactory *funFact;
  //AIDA::IFitFactory      *fitFact;
  //AIDA::IFunction        *exponFun;
  //AIDA::IFunction        *gaussFun;
  //AIDA::IFitter          *e_fitter;
  ///AIDA::IFitResult       *fitResult;
};
#endif
#endif



