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
#include "Par03Hit.hh"

#include "G4VisAttributes.hh"
#include "G4Tubs.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VVisManager.hh"
#include "G4LogicalVolume.hh"

G4ThreadLocal G4Allocator<Par03Hit>* Par03HitAllocator;

Par03Hit::Par03Hit()
  : G4VHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par03Hit::~Par03Hit() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par03Hit::Par03Hit(const Par03Hit& aRight)
  : G4VHit()
{
  fEdep   = aRight.fEdep;
  fZId    = aRight.fZId;
  fRhoId  = aRight.fRhoId;
  fPhiId  = aRight.fPhiId;
  fTime   = aRight.fTime;
  fPos    = aRight.fPos;
  fRot    = aRight.fRot;
  fType   = aRight.fType;
  fLogVol = aRight.fLogVol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const Par03Hit& Par03Hit::operator=(const Par03Hit& aRight)
{
  fEdep   = aRight.fEdep;
  fZId    = aRight.fZId;
  fRhoId  = aRight.fRhoId;
  fPhiId  = aRight.fPhiId;
  fTime   = aRight.fTime;
  fPos    = aRight.fPos;
  fRot    = aRight.fRot;
  fType   = aRight.fType;
  fLogVol = aRight.fLogVol;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int Par03Hit::operator==(const Par03Hit& aRight) const
{
  return (fRhoId == aRight.fRhoId && fPhiId == aRight.fPhiId &&
          fZId == aRight.fZId);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03Hit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  // Hits can be filtered out in visualisation
  if(!pVVisManager->FilterHit(*this))
    return;
  // Do not draw empty hits
  if(fEdep < 0)
    return;
  if(pVVisManager)
  {
    G4Transform3D trans(fRot, fPos);
    G4VisAttributes attribs;
    // Create default dimensions
    G4Tubs solid("draw", 0, 1 * cm, 1 * cm, 0, 0.05 * CLHEP::pi);
    if(fLogVol)
    {
      const G4VisAttributes* pVA = fLogVol->GetVisAttributes();
      if(pVA)
        attribs = *pVA;
      // Cannot use directly fLogVol due to rho parametrisation (change of
      // solid!) Recalculation of radius is needed
      solid     = *dynamic_cast<G4Tubs*>(fLogVol->GetSolid());
      double dR = solid.GetOuterRadius() - solid.GetInnerRadius();
      solid.SetInnerRadius(solid.GetInnerRadius() + fRhoId * dR);
      solid.SetOuterRadius(solid.GetOuterRadius() + fRhoId * dR);
    }
    // Set colours depending on type of hit (full/fast sim)
    G4double colR = fType == 0 ? 0 : 1;
    G4double colG = fType == 0 ? 1 : 0;
    G4double colB = 0;
    G4Colour colour(colR, colG, colB, 0.5);
    attribs.SetColour(colour);
    attribs.SetForceSolid(true);
    pVVisManager->Draw(solid, attribs, trans);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::map<G4String, G4AttDef>* Par03Hit::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String, G4AttDef>* store =
    G4AttDefStore::GetInstance("Par03Hit", isNew);
  if(isNew)
  {
    (*store)["HitType"] =
      G4AttDef("HitType", "Hit Type", "Physics", "", "G4String");
    (*store)["Energy"] = G4AttDef("Energy", "Energy Deposited", "Physics",
                                  "G4BestUnit", "G4double");
    (*store)["Time"] =
      G4AttDef("Time", "Time", "Physics", "G4BestUnit", "G4double");
    (*store)["Pos"] =
      G4AttDef("Pos", "Position", "Physics", "G4BestUnit", "G4ThreeVector");
  }
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4AttValue>* Par03Hit::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;
  values->push_back(G4AttValue("HitType", "HadPar03Hit", ""));
  values->push_back(G4AttValue("Energy", G4BestUnit(fEdep, "Energy"), ""));
  values->push_back(G4AttValue("Time", G4BestUnit(fTime, "Time"), ""));
  values->push_back(G4AttValue("Pos", G4BestUnit(fPos, "Length"), ""));
  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03Hit::Print()
{
  std::cout << "\tHit " << fEdep / MeV << " MeV at " << fPos / cm
            << " cm (R,phi,z)= (" << fRhoId << ", " << fPhiId << ", " << fZId
            << "), " << fTime << " ns" << std::endl;
}
