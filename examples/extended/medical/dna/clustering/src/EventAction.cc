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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"

#include "Analysis.hh"
#include "ClusteringAlgo.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction():G4UserEventAction()
{
  //default parameter values
  fEdep=0.;

  // Create clustering algorithm
  // These default values have been tuned for the Physics List G4EmDNAPhysics
  // to reproduce data published by:
  // Francis et al. 2011 Comput. Meth. Programs. Biomed. 2011 101(3)
  fpClustering = new ClusteringAlgo(3.3*nanometer,2,0.2,5*eV,37.5*eV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  delete fpClustering;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction( const G4Event*)
{
  fEdep=0.;
  fpClustering->Purge();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction( const G4Event*)
{
  std::map<G4int,G4int> sizeDistribution = fpClustering->RunClustering();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  analysisManager->FillH1(1,fpClustering->GetSSB());
  analysisManager->FillH1(2,fpClustering->GetComplexSSB());
  analysisManager->FillH1(3,fpClustering->GetDSB());

  while ( !sizeDistribution.empty() )
  {
    analysisManager->FillH1(4,
                            sizeDistribution.begin()->first,
                            sizeDistribution.begin()->second);
    sizeDistribution.erase(sizeDistribution.begin());
  } 

  analysisManager->FillH1(5,
                          (fEdep/joule)/
                          (G4LogicalVolumeStore::GetInstance()->
                              GetVolume("Target")->GetMass()/kg)
  );
}
