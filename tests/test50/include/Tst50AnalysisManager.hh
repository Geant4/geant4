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
//    **********************************
//    *                                *
//    *      Tst50AnalysisManager.hh  *
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
}


class Tst50AnalysisManager { 

public:
  
  ~Tst50AnalysisManager();
  static Tst50AnalysisManager* getInstance();
  void book(G4String);
  void bookHistograms(); 
  void AttenuationGammaCoeffiecient(G4int,G4double,G4double,G4double);
  void StoppingPower(G4int,G4double,G4double);
  void CSDARange(G4int,G4double,G4double);
  void ParticleTransmission(G4int,G4double,G4double,G4double,G4double,G4double);  void TransmittedEnergy(G4int,G4double,G4double); 
  void FillEnergyDeposit(G4double);
  void finish();

private:

  static Tst50AnalysisManager* instance;
  Tst50AnalysisManager();

  AIDA::IAnalysisFactory*  aFact; 
  AIDA::ITreeFactory*      treeFact;
  AIDA::ITree*             theTree;

  AIDA::IDataPointSetFactory *  dataPointFactory; 
  AIDA::IHistogramFactory*     histogramFactory;

  AIDA::IDataPointSet *  stoppingPowerDataPoint;  
  AIDA::IDataPointSet *  CSDARangeDataPoint;
  AIDA::IDataPointSet *  particleTransmissionDataPoint;
  AIDA::IDataPointSet * particleTransmissionEnergyDataPoint;
  AIDA::IDataPointSet *  gammaAttenuationCoefficientDataPoint;
  AIDA::IHistogram1D* histogramEnergyDeposit;
};
#endif
#endif




