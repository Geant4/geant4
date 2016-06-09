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
// ----------------------------------------------------------------------------
// $Id: HadrontherapyAnalysisManager.hh; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova, Genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#ifdef G4ANALYSIS_USE
#ifndef HADRONTHERAPYANALYSISMANAGER_HH
#define HADRONTHERAPYANALYSISMANAGER_HH 1

#include "globals.hh"
# include <AIDA/AIDA.h>

namespace AIDA{
  class ITree; 
  class IAnalysisFactory;
  class ITreeFactory;
};

class HadrontherapyAnalysisManager
{
private:
  HadrontherapyAnalysisManager();

public:
  ~HadrontherapyAnalysisManager();
  
  static HadrontherapyAnalysisManager* getInstance();
  
  void book();
  // Book the histograms and ntuples in a .hbk file
  
  void FillEnergyDeposit(G4int voxelXId, G4int voxelYId, G4int voxelZId, 
                         G4double energyDeposit);
  // Fill the ntuple with the energy deposit in the phantom 

  void BraggPeak(G4int, G4double);
  // Fill 1D histogram with the Bragg peak in the phantom
  
  void finish();
  // Close the .hbk file with the histograms and the ntuples

private:
  static HadrontherapyAnalysisManager* instance;
  AIDA::IAnalysisFactory* aFact;
  AIDA::ITree* theTree; 
  AIDA::IHistogramFactory *histFact;
  AIDA::ITupleFactory *tupFact;
  AIDA::IHistogram1D *h1;
  AIDA::ITuple *ntuple;
};
#endif
#endif



