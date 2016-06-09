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
// $Id: RemSimAnalysisManager.hh,v 1.7 2004/11/23 11:43:21 guatelli Exp $
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
 
#ifdef G4ANALYSIS_USE
#ifndef REMSIMANALYSISMANAGER_HH
#define REMSIMANALYSISMANAGER_HH 

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
# include <AIDA/AIDA.h>
#include "RemSimAnalysisMessenger.hh"

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

class RemSimAnalysisMessenger;
class RemSimAnalysisManager { 

public:
  
  ~RemSimAnalysisManager();
  static RemSimAnalysisManager* getInstance();
  void book(); // booking the hbook file

  void energyDepositStore(G4int, G4double); 
  // Collect the energy deposit in the phantom 
                                           
  void primaryParticleEnergyDistribution(G4double);
  // Initail energy of primary particles

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

  void particleShape(G4double, G4double);

  void energyDepShape(G4double, G4double, G4double);

  void SetFormat(G4String);

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
  AIDA::IHistogram2D* shape;
  AIDA::IHistogram2D* energyShape;
  
  RemSimAnalysisMessenger* messenger;
  G4String fileFormat; 
};
#endif
#endif




