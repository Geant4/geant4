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
/// \file CalorimeterSD.cc
/// \brief Implementation of the CaTS::CalorimeterSD class

// Geant4 headers
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
// project headers
#include "CalorimeterSD.hh"
#include "ConfigurationManager.hh"

CalorimeterSD::CalorimeterSD(G4String name)
  : G4VSensitiveDetector(name)
{
  G4String HCname = name + "_HC";
  collectionName.insert(HCname);
  G4cout << collectionName.size() << "   CalorimeterSD name:  " << name
         << " collection Name: " << HCname << G4endl;
  fHCID   = -1;
  verbose = ConfigurationManager::getInstance()->isEnable_verbose();
}

void CalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  fCalorimeterHitsCollection =
    new CalorimeterHitsCollection(SensitiveDetectorName, collectionName[0]);
  if(fHCID < 0)
  {
    if(verbose)
      G4cout << "CalorimeterSD::Initialize:  " << SensitiveDetectorName << "   "
             << collectionName[0] << G4endl;
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  hce->AddHitsCollection(fHCID, fCalorimeterHitsCollection);
}

G4bool CalorimeterSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  // auto can be used for all results from retrieving functions
  auto edep = aStep->GetTotalEnergyDeposit() / CLHEP::MeV;
  if(edep == 0.)
    return false;
  auto time       = aStep->GetPreStepPoint()->GetGlobalTime() / CLHEP::ns;
  auto touch      = aStep->GetPreStepPoint()->GetTouchable();
  auto cellpos    = touch->GetTranslation();
  unsigned int ID = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();
  auto theTrack   = aStep->GetTrack();
  auto particleType = theTrack->GetDefinition()->GetParticleName();
  //
  //  check if this cell has been hit before
  // fCalorimeterHitsCollection
  for(unsigned int j = 0; j < fCalorimeterHitsCollection->entries(); j++)
  {
    CalorimeterHit* aPreviousHit = (*fCalorimeterHitsCollection)[j];
    if(ID == aPreviousHit->GetId())
    {
      aPreviousHit->SetEdep(aStep->GetTotalEnergyDeposit() +
                            aPreviousHit->GetEdep());
      if((particleType == "e+") || (particleType == "gamma") ||
         (particleType == "e-"))
      {
        aPreviousHit->Setem_Edep(edep + aPreviousHit->Getem_Edep());
      }
      return true;
    }
  }
  //
  // otherwise create a new hit:
  //
  CalorimeterHit* newHit;
  if((particleType == "e+") || (particleType == "gamma") ||
     (particleType == "e-"))
  {
    newHit = new CalorimeterHit(ID, edep, edep, time, cellpos);
  }
  else
  {
    newHit = new CalorimeterHit(ID, edep, 0.0, time, cellpos);
  }
  fCalorimeterHitsCollection->insert(newHit);
  return true;
}

void CalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  G4int NbHits = fCalorimeterHitsCollection->entries();
  if(verbose)
    G4cout << " Number of CalorimeterHits:  " << NbHits << G4endl;
}
