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
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// ScintSD (scintillator sensitive detector definition) program
// --------------------------------------------------------------

#include "DMXScintSD.hh"

#include "DMXScintHit.hh"
#include "DMXDetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Ions.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXScintSD::DMXScintSD(G4String name,DMXDetectorConstruction* DMXSD) 
  :G4VSensitiveDetector(name),DMXDetector(DMXSD)
{
  G4String HCname="scintillatorCollection";
  collectionName.insert(HCname);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXScintSD::~DMXScintSD(){ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXScintSD::Initialize(G4HCofThisEvent* HCE)
{
  scintillatorCollection = new DMXScintHitsCollection
    (SensitiveDetectorName,collectionName[0]);

  HitID = -1;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool DMXScintSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4ParticleDefinition* particleType = aStep->GetTrack()->GetDefinition();
  G4String particleName = particleType->GetParticleName();

  G4double stepl = 0.;
  if (particleType->GetPDGCharge() != 0.)
    stepl = aStep->GetStepLength();
  
  if ((edep==0.)&&(stepl==0.)) return false;      


  // fill in hit
  DMXScintHit* newHit = new DMXScintHit();
  newHit->SetEdep(edep);
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  newHit->SetTime(aStep->GetPreStepPoint()->GetGlobalTime());
  newHit->SetParticle(particleName);
  newHit->SetParticleEnergy(aStep->GetPreStepPoint()->GetKineticEnergy() );

  HitID = scintillatorCollection->insert(newHit);
  
  return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXScintSD::EndOfEvent(G4HCofThisEvent* HCE)
{

  G4String HCname = collectionName[0];
  static G4int HCID = -1;
  if(HCID<0)
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(HCname);
  HCE->AddHitsCollection(HCID,scintillatorCollection);

  G4int nHits = scintillatorCollection->entries();
  if (verboseLevel>=1)
    G4cout << "     LXe collection: " <<  nHits << " hits" << G4endl;
  if (verboseLevel>=2)
    scintillatorCollection->PrintAllHits();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXScintSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXScintSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXScintSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

