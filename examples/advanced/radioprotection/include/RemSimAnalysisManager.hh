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
  void book();
  void energyDeposit1(G4double,G4double);
  void energyDeposit2(G4double,G4double);
  void energyDeposit3(G4double,G4double);
  void energyDepositStore(G4int, G4double);
  void leptonsEnergySpectrum1(G4double);
  void hadronEnergySpectrum1(G4double);
  void gammaEnergySpectrum1(G4double);
  void leptonsEnergySpectrum2(G4double);
  void hadronEnergySpectrum2(G4double);
  void gammaEnergySpectrum2(G4double);
  void leptonsEnergySpectrum3(G4double);
  void hadronEnergySpectrum3(G4double);
  void gammaEnergySpectrum3(G4double);
  void neutronEnergyDistribution(G4double);
  void photonEnergyDistribution(G4double);
  void electronEnergyDistribution(G4double);
  void hadronEnergyDistribution(G4double);
  void primaryParticleEnergyDistribution(G4double);
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
  AIDA::IHistogram2D* histo1;
  AIDA::IHistogram2D* histo2;
  AIDA::IHistogram2D* histo3;
  AIDA::IHistogram1D* trasmission1;  
  AIDA::IHistogram1D* trasmission2; 
  AIDA::IHistogram1D* trasmission3;  
  AIDA::IHistogram1D* trasmission12;  
  AIDA::IHistogram1D* trasmission22; 
  AIDA::IHistogram1D* trasmission32; 
  AIDA::IHistogram1D* trasmission13;  
  AIDA::IHistogram1D* trasmission23; 
  AIDA::IHistogram1D* trasmission33;
  AIDA::IHistogram1D* neutron;
  AIDA::IHistogram1D* photon;
  AIDA::IHistogram1D* electron; 
  AIDA::IHistogram1D* hadron;
  AIDA::IHistogram1D* primary; 

};
#endif
#endif




