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
// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file DRCalorimeterSD.cc
/// \brief Implementation of the CaTS::DRCalorimeterSD class

// Geant4 headers
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4EventManager.hh"
#include "G4Cerenkov.hh"
// project headers
#include "DRCalorimeterSD.hh"
#include "ConfigurationManager.hh"

DRCalorimeterSD::DRCalorimeterSD(G4String name)
  : G4VSensitiveDetector(name)
  , fDRCalorimeterHitsCollection(0)
  , fHCID(0)
{
  G4String HCname = name + "_HC";
  collectionName.insert(HCname);
  G4cout << collectionName.size() << "   DRCalorimeterSD name:  " << name
         << " collection Name: " << HCname << G4endl;
  fHCID   = -1;
  verbose = ConfigurationManager::getInstance()->isEnable_verbose();
}

DRCalorimeterSD::~DRCalorimeterSD() {}

void DRCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  fDRCalorimeterHitsCollection =
    new DRCalorimeterHitsCollection(SensitiveDetectorName, collectionName[0]);
  if(fHCID < 0)
  {
    if(verbose)
      G4cout << "DRCalorimeterSD::Initialize:  " << SensitiveDetectorName
             << "   " << collectionName[0] << G4endl;
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  hce->AddHitsCollection(fHCID, fDRCalorimeterHitsCollection);
}

G4bool DRCalorimeterSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit() / CLHEP::MeV;
  if(edep == 0.)
    return false;
  G4double time = aStep->GetPreStepPoint()->GetGlobalTime() / CLHEP::ns;
  const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();
  G4ThreeVector cellpos     = touch->GetTranslation();
  unsigned int ID = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();
  G4Track* theTrack          = aStep->GetTrack();
  G4String particleType      = theTrack->GetDefinition()->GetParticleName();
  unsigned int NCerenPhotons = 0;
  G4SteppingManager* fpSteppingManager = G4EventManager::GetEventManager()
                                           ->GetTrackingManager()
                                           ->GetSteppingManager();
  G4StepStatus stepStatus = fpSteppingManager->GetfStepStatus();
  if(stepStatus != fAtRestDoItProc)
  {
    G4ProcessVector* procPost = fpSteppingManager->GetfPostStepDoItVector();
    size_t MAXofPostStepLoops = fpSteppingManager->GetMAXofPostStepLoops();
    for(size_t i3 = 0; i3 < MAXofPostStepLoops; i3++)
    {
      if((*procPost)[i3]->GetProcessName() == "Cerenkov")
      {
        G4Cerenkov* proc = (G4Cerenkov*) (*procPost)[i3];
        NCerenPhotons += proc->GetNumPhotons();
      }
    }
  }
  //
  //  check if this cell has been hit before
  // fDRCalorimeterHitsCollection
  for(unsigned int j = 0; j < fDRCalorimeterHitsCollection->entries(); j++)
  {
    DRCalorimeterHit* aPreviousHit = (*fDRCalorimeterHitsCollection)[j];
    if(ID == aPreviousHit->GetId())
    {
      aPreviousHit->SetEdep(aStep->GetTotalEnergyDeposit() +
                            aPreviousHit->GetEdep());
      aPreviousHit->SetNceren(aPreviousHit->GetNceren() + NCerenPhotons);
      if((particleType == "e+") || (particleType == "gamma") ||
         (particleType == "e-"))
      {
        aPreviousHit->SetEm_Edep(edep + aPreviousHit->GetEm_Edep());
      }
      return true;
    }
  }
  //
  // otherwise create a new hit:
  //
  DRCalorimeterHit* newHit;
  if((particleType == "e+") || (particleType == "gamma") ||
     (particleType == "e-"))
  {
    newHit = new DRCalorimeterHit(ID, edep, edep, NCerenPhotons, time, cellpos);
  }
  else
  {
    newHit = new DRCalorimeterHit(ID, edep, 0.0, time, NCerenPhotons, cellpos);
  }
  fDRCalorimeterHitsCollection->insert(newHit);
  return true;
}

void DRCalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  G4int NbHits = fDRCalorimeterHitsCollection->entries();
  if(verbose)
    G4cout << " Number of DRCalorimeterHits:  " << NbHits << G4endl;
}
