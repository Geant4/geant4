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
//    **********************************
//    *                                *
//    *      RemSimAnalysisManager.hh  *
//    *                                *
//    **********************************
// 
//
// the class Analysis creates and managed histograms and ntuples
///
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
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
};

class RemSimAnalysisManager { 

public:
  
  ~RemSimAnalysisManager();
  static RemSimAnalysisManager* getInstance();
  void book(); // booking the hbook file

  void energyDepositStore(G4int, G4double); 
  // Collect the energy deposit in the phantom 
                                           
  void primaryParticleEnergyDistribution(G4double);
  // Energy of primary particles

  void SecondaryEnergyDeposit(G4int, G4double);
  // Energy deposit given by secondary particles in the phantom

  void PrimaryInitialEnergyIn(G4double);
  // Initial energy of primary particles impinging on the phantom

  void PrimaryInitialEnergyOut(G4double);
  // Initial energy of primary particles outgoing the phantom 

  void PrimaryEnergyIn(G4double);
  // Energy of primary particles impinging on the phantom

  void PrimaryEnergyOut(G4double);
  // Energy of primary particles outgoing the phantom 

  void finish();

private:

  static RemSimAnalysisManager* instance;
  RemSimAnalysisManager();

  AIDA::IAnalysisFactory*  aFact; 
  AIDA::ITreeFactory*      treeFact;
  AIDA::ITree*             theTree;

  AIDA::IDataPointSetFactory *  dataPointFactory; 
  AIDA::IHistogramFactory*     histogramFactory;

  AIDA::IDataPointSet *  dataPoint; 
  AIDA::IHistogram1D* energyDeposit; 
  AIDA::IHistogram1D* primary;  
  AIDA::IHistogram1D* secondaryDeposit;
  AIDA::IHistogram1D* primaryInitialE;
  AIDA::IHistogram1D* primaryInitialEout;
  AIDA::IHistogram1D* initialE;
  AIDA::IHistogram1D* initialEout;
};
#endif
#endif




