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
//
// $Id: DMXAnalysisManager.hh
// GEANT4 tag $Name:
//
// Author: Alex Howard (a.s.howard@ic.ac.uk)
//
// History:
// -----------
//  16 Jan 2002  Alex Howard   Created
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

#include "AIDA/IPlotterFactory.h"
#include "AIDA/IPlotter.h"

#include "AIDA/ITupleFactory.h"
#include "AIDA/ITuple.h"

#include "AIDA/IManagedObject.h"


class IAnalysisFactory;
class ITree;
class IHistogramFactory;
class ITupleFactory;
class ITuple;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class DMXAnalysisManager
{
public:
 
  virtual ~DMXAnalysisManager();
  
  //  void book();
  void book(G4String);
  
  void finish();
  
  //fill histograms with SHC (Scint Hits) data from DMXEventAction
  void analyseScintHits(G4int event_id, G4double energy_pri, G4double totEnergy, G4int S_hits, G4double firstLXeHitTime, G4int P_hits, G4double aveTimePmtHits, G4String firstparticleName, G4double firstParticleE, G4bool gamma_ev, G4bool neutron_ev, G4bool positron_ev, G4bool electron_ev, G4bool other_ev, long seed1, long seed2);

  //fill histograms with PHC (Pmt Hits) data from DMXEventAction
  void analysePMTHits(G4int, G4int, G4double, G4double, G4double);
  
  //fill histograms with data from DMXParticleSource / secondary history
  void analyseParticleSource(G4double, G4String);

  //fill histograms with data from DMXPrimaryGenerator
  void analysePrimaryGenerator(G4double);
  
  //fill histograms with data from DMXPrimaryGenerator
  void HistFirstTime(G4double);
  
  //Method to interactively display histograms
  void PlotHistos();
  
  //method to call to create an instance of this class
  static DMXAnalysisManager* getInstance();
 

private:
  
  //private constructor in order to create a singleton
  DMXAnalysisManager();
 
  static DMXAnalysisManager* instance;
  
  // Quantities for the ntuple


  G4int event_id; G4double energy_pri; G4double totEnergy; G4int S_hits; G4double firstLXeHitTime; G4int P_hits; G4double aveTimePmtHits; G4String firstparticleName; G4double firstParticleE; G4bool gamma_ev; G4bool neutron_ev; G4bool positron_ev;G4bool electron_ev;G4bool other_ev; long seed1; long seed2;


  G4int i; G4double x; G4double y; G4double z;

  G4double energy; G4double name;

  G4double time;


  IAnalysisFactory  *af;
  ITree             *tree;
  IHistogramFactory *hf;
  ITupleFactory     *tpf;
  IPlotterFactory   *pf;
  IPlotter          *plotter;

};
#endif
#endif



