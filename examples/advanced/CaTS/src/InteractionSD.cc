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
/// \file InteractionSD.cc
/// \brief Implementation of the CaTS::InteractionSD class

// Geant4 headers:
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
// project headers
#include "InteractionSD.hh"
#include "ParticleChange.hh"

InteractionSD::InteractionSD(G4String name)
  : G4VSensitiveDetector(name)
{
  G4String HCname = name + "_HC";
  collectionName.insert(HCname);
  G4cout << collectionName.size() << "   InteractionSD name:  " << name
         << " collection Name: " << HCname << G4endl;
  fHCID = -1;
}

InteractionSD::~InteractionSD()
{
  delete fFirstInter;
  delete fOtherInter;
}

void InteractionSD::Initialize(G4HCofThisEvent* HCE)
{
  G4cout << "Hits Collection capacity:  " << HCE->GetCapacity() << G4endl;
  fInteractionHitsCollection =
    new InteractionHitsCollection(SensitiveDetectorName, collectionName[0]);
  if(fHCID < 0)
  {
    G4cout << "InteractionSD::Initialize:  " << SensitiveDetectorName << "   "
           << collectionName[0] << G4endl;
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection(fHCID, fInteractionHitsCollection);
  fFirstInter = new ParticleChange(true);
  fOtherInter = new ParticleChange();
}

G4bool InteractionSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4int nsc = fFirstInter->GetNumberOfSecondaries();
  if(nsc > 0)
  {
    for(G4int i = 0; i < nsc; i++)
    {
      delete fFirstInter->GetSecondary(i);
    }
    fFirstInter->Clear();
  }
  nsc = fOtherInter->GetNumberOfSecondaries();
  if(nsc > 0)
  {
    for(G4int i = 0; i < nsc; i++)
    {
      delete fOtherInter->GetSecondary(i);
    }
    fOtherInter->Clear();
  }
  const std::vector<const G4Track*>* secs = aStep->GetSecondaryInCurrentStep();
  G4int nsec                              = secs->size();
  for(G4int i = 0; i < nsec; i++)
  {
    G4Track* tr = new G4Track(*((*secs)[i]));
    if(aStep->GetTrack()->GetTrackStatus() != fAlive)  // track looses identity
    {
      if(aStep->GetTrack()->GetParentID() == 0)  // primary track
      {
        fFirstInter->AddSecondary(tr);
      }
      else  // secondary track, and it's also looses identity (re-interaction)
      {
        fOtherInter->AddSecondary(tr);
      }
    }
  }  // end loop over secondaries
  G4int NSec = fFirstInter->GetNumberOfSecondaries();
  if(NSec > 0)
  {
    const G4DynamicParticle* sec = 0;
    for(G4int i = 0; i < NSec; i++)
    {
      sec = fFirstInter->GetSecondary(i)->GetDynamicParticle();
      const G4String& pname  = sec->GetDefinition()->GetParticleName();
      G4double pmom          = (sec->GetTotalMomentum()) / GeV;
      G4double Ekin          = (sec->GetKineticEnergy()) / GeV;
      G4double theta         = (sec->GetMomentum()).theta();
      InteractionHit* newHit = new InteractionHit(pname, pmom, Ekin, theta);
      fInteractionHitsCollection->insert(newHit);
    }
  }
  return true;
}
