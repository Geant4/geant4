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
//    *      RemSimAnalysisManager.hh  *
//    *                                *
//    **********************************
//
// $Id: RemSimAnalysisManager.hh,v 1.12 2009-04-08 13:25:37 cirrone Exp $
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

class RemSimAnalysisMessenger;
class RemSimAnalysisManager { 

public:
  
  ~RemSimAnalysisManager();
  static RemSimAnalysisManager* getInstance();
  void book(); // booking the hbook file

  void energyDepositStore(G4int, G4double); 
  // Collect the energy deposit in the phantom 
  // deriving from primary and secondary particles                                         

  void primaryParticleEnergyDistribution(G4double);
  // Initial energy of primary particles

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
  // Position of the energy deposits in the phantom, projected in the plane perpendicular
  // to the direction of the primary particles (z axis)

  void energyDepShape(G4double, G4double, G4double);
  // Energy deposits in the phantom, projected in the plane perpendicular
  // to the direction of the primary particles (z axis)

  void SetFormat(G4String);
  // At the moment it is not possible to switch the formta of the file
  // Number of secondary particles produced in the phantom
  void SecondaryInPhantom(G4int); 

  // Energy spectrum of secondary particles originated in the phantom
  void SecondaryProtonInPhantom(G4double);
  void SecondaryNeutronInPhantom(G4double);
  void SecondaryPionInPhantom(G4double);
  void SecondaryAlphaInPhantom(G4double);
  void SecondaryPositronInPhantom(G4double);
  void SecondaryElectronInPhantom(G4double);
  void SecondaryGammaInPhantom(G4double);
  void SecondaryMuonInPhantom(G4double);
  void SecondaryOtherInPhantom(G4double);
  void SecondaryNeutrinoInPhantom(G4double);
  void SecondaryProtonInPhantomSlice(G4double);
  void SecondaryNeutronInPhantomSlice(G4double);
  void SecondaryPionInPhantomSlice(G4double);
  void SecondaryAlphaInPhantomSlice(G4double);
  void SecondaryPositronInPhantomSlice(G4double);
  void SecondaryElectronInPhantomSlice(G4double);
  void SecondaryGammaInPhantomSlice(G4double);
  void SecondaryMuonInPhantomSlice(G4double);
  void SecondaryOtherInPhantomSlice(G4double);

  // Number of secondary particles reaching the phantom
  void SecondaryReachingThePhantom(G4int);

  // Energy spectrum of secondary particles reaching the phantom
  void SecondaryProtonReachingThePhantom(G4double);
  void SecondaryNeutronReachingThePhantom(G4double);
  void SecondaryPionReachingThePhantom(G4double);
  void SecondaryAlphaReachingThePhantom(G4double);
  void SecondaryPositronReachingThePhantom(G4double);
  void SecondaryElectronReachingThePhantom(G4double);
  void SecondaryGammaReachingThePhantom(G4double);
  void SecondaryMuonReachingThePhantom(G4double);
  void SecondaryOtherReachingThePhantom(G4double);
  
 // Number of secondary particles produced in the vehicle
  void SecondaryInVehicle(G4int);

  void finish();

private:
  static RemSimAnalysisManager* instance;
  RemSimAnalysisManager();

  AIDA::IAnalysisFactory*  aFact; 
  AIDA::ITreeFactory*      treeFact;
  AIDA::ITree*             theTree;

  AIDA::IDataPointSetFactory * dataPointFactory; 
  AIDA::IHistogramFactory*     histogramFactory;

  AIDA::IDataPointSet* dataPoint; 
  AIDA::IHistogram1D* energyDeposit; 
  AIDA::IHistogram1D* primary;  
  AIDA::IHistogram1D* secondaryDeposit; 
  AIDA::IHistogram1D* primaryInitialE;
  AIDA::IHistogram1D* primaryInitialEout;
  AIDA::IHistogram1D* initialE;
  AIDA::IHistogram1D* initialEout;
  AIDA::IHistogram2D* shape;
  AIDA::IHistogram2D* energyShape;

  AIDA::IHistogram1D* histo_secondary_phantom;
  AIDA::IHistogram1D* histo_secondary;
  AIDA::IHistogram1D* histo_proton;
  AIDA::IHistogram1D* histo_neutron;
  AIDA::IHistogram1D* histo_pion;
  AIDA::IHistogram1D* histo_alpha;
  AIDA::IHistogram1D* histo_positron;
  AIDA::IHistogram1D* histo_electron;
  AIDA::IHistogram1D* histo_gamma;
  AIDA::IHistogram1D* histo_muon;
  AIDA::IHistogram1D* histo_other;
  AIDA::IHistogram1D* histo_neutrino;
  AIDA::IHistogram1D* histo_proton_reaching;
  AIDA::IHistogram1D* histo_neutron_reaching;
  AIDA::IHistogram1D* histo_pion_reaching;
  AIDA::IHistogram1D* histo_alpha_reaching;
  AIDA::IHistogram1D* histo_positron_reaching;
  AIDA::IHistogram1D* histo_electron_reaching;
  AIDA::IHistogram1D* histo_gamma_reaching;
  AIDA::IHistogram1D* histo_muon_reaching;
  AIDA::IHistogram1D* histo_other_reaching;

  AIDA::IHistogram1D* histo_proton_slice;
  AIDA::IHistogram1D* histo_neutron_slice;
  AIDA::IHistogram1D* histo_pion_slice;
  AIDA::IHistogram1D* histo_alpha_slice; 
  AIDA::IHistogram1D* histo_positron_slice;
  AIDA::IHistogram1D* histo_electron_slice;
  AIDA::IHistogram1D* histo_gamma_slice;
  AIDA::IHistogram1D* histo_muon_slice;
  AIDA::IHistogram1D* histo_other_slice;
  AIDA::IHistogram1D* histo_vehicle;
 
  G4String fileFormat; 
};
#endif
#endif




