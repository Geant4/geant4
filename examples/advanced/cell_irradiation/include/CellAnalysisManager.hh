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
//    **************************************
//    *                                    *
//    *       CellAnalysisManager.hh       *
//    *                                    *
//    **************************************
// 

// the class Analysis creates and managed histograms and ntuples
///
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//	   Barbara Mascialino (Barbara.Mascialino@ge.infn.it)
//
// History:
// -----------
// 20 September  2006   S. Guatelli, B. Mascialino   1st implementation
//
// -------------------------------------------------------------------
 
#ifdef G4ANALYSIS_USE
#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
# include <AIDA/AIDA.h>

namespace AIDA 
{
  class ITree;
  class IHistogramFactory;
  class IAnalysisFactory;
  class IDataPoint;
  class ITupleFactory;
  class ITuple;
  class ITreeFactory;
}


class CellAnalysisManager { 

public:
  
  ~CellAnalysisManager();
  static CellAnalysisManager* getInstance();

  // Creating the hbk file at the beginning of the simulation
  void book(G4String);
  
  // Fill with the initial energy of primary particles
  void primaryparticle_energy(G4double);
  
  // Fill with the energy deposit (MeV) in the plane, perpendicular 
  // to z
  void FillEnergyDeposit(G4double, G4double, G4double); 

  // Fiil with the energy the primary particles have when they outgo
  // the target
  void primary_energy_outgoing(G4double);

  void FillProfile(G4double, G4double);

  // Closes the hbk file with the results stored at the end of the simulation
  void finish();

private:

  static CellAnalysisManager* instance;
  CellAnalysisManager();

  AIDA::IAnalysisFactory*  aFact; 
  AIDA::ITreeFactory*      treeFact;
  AIDA::ITree*             theTree;
  AIDA::IHistogramFactory*     histogramFactory;

  AIDA::IHistogram1D* primary_energy;
  AIDA::IHistogram2D* histogramEnergyDeposit;
  AIDA::IHistogram1D* primary_outgoing_energy;
  AIDA::IHistogram1D* profile;
};
#endif
#endif




