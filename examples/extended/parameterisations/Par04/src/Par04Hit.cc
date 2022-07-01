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
#include "Par04Hit.hh"
#include <CLHEP/Units/SystemOfUnits.h>     // for mm, pi, MeV, cm, rad
#include <CLHEP/Vector/ThreeVector.h>      // for operator/, operator<<, Hep...
#include <G4RotationMatrix.hh>             // for G4RotationMatrix
#include <G4String.hh>                     // for G4String
#include <G4ThreeVector.hh>                // for G4ThreeVector
#include <G4Transform3D.hh>                // for G4Transform3D
#include <G4VHit.hh>                       // for G4VHit
#include <algorithm>                       // for max
#include <iostream>                        // for operator<<, basic_ostream:...
#include <string>                          // for operator<
#include "G4AttDef.hh"                     // for G4AttDef
#include "G4AttDefStore.hh"                // for GetInstance
#include "G4AttValue.hh"                   // for G4AttValue
#include "G4Colour.hh"                     // for G4Colour
#include "G4SystemOfUnits.hh"              // for mm, MeV, cm, rad
#include "G4Tubs.hh"                       // for G4Tubs
#include "G4UnitsTable.hh"                 // for G4BestUnit
#include "G4VVisManager.hh"                // for G4VVisManager
#include "G4VisAttributes.hh"              // for G4VisAttributes
#include <cmath>                           // for log10
template <class Type> class G4Allocator;

G4ThreadLocal G4Allocator<Par04Hit>* Par04HitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04Hit::Par04Hit()
  : G4VHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04Hit::~Par04Hit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04Hit::Par04Hit(const Par04Hit& aRight)
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

const Par04Hit& Par04Hit::operator=(const Par04Hit& aRight)
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

int Par04Hit::operator==(const Par04Hit& aRight) const
{
  return (fRhoId == aRight.fRhoId && fPhiId == aRight.fPhiId && fZId == aRight.fZId);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04Hit::Draw()
{
  /// Arbitrary size corresponds to the example macros
  G4ThreeVector meshSize(2.325 * mm, 2 * CLHEP::pi / 50. * CLHEP::rad, 3.4 * mm);
  G4int numPhiCells           = CLHEP::pi * 2. / meshSize.y();
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  // Hits can be filtered out in visualisation
  if(!pVVisManager->FilterHit(*this))
    return;
  // Do not draw empty hits
  if(fEdep <= 0)
    return;
  // Do not plot if default values were not changed
  if(fRhoId == -1 && fZId == -1 && fPhiId == -1)
    return;
  if(pVVisManager)
  {
    G4Transform3D trans(fRot, fPos);
    G4VisAttributes attribs;
    G4Tubs solid("draw", fRhoId * meshSize.x(), (fRhoId + 1) * meshSize.x(), meshSize.z() / 2.,
                 (-numPhiCells / 2. + fPhiId) * meshSize.y(), meshSize.y());
    // Set colours depending on type of hit (full/fast sim)
    G4double colR = fType == 0 ? 0 : 1;
    G4double colG = fType == 0 ? 1 : 0;
    G4double colB = 0;
    // Set transparency depending on the energy
    // Arbitrary formula
    G4double alpha = 2 * std::log10(fEdep + 1);
    G4cout << "alpha = " << alpha << G4endl;
    G4Colour colour(colR, colG, colB, alpha);
    attribs.SetColour(colour);
    attribs.SetForceSolid(true);
    pVVisManager->Draw(solid, attribs, trans);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::map<G4String, G4AttDef>* Par04Hit::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String, G4AttDef>* store = G4AttDefStore::GetInstance("Par04Hit", isNew);
  if(isNew)
  {
    (*store)["HitType"] = G4AttDef("HitType", "Hit Type", "Physics", "", "G4String");
    (*store)["Energy"] =
      G4AttDef("Energy", "Energy Deposited", "Physics", "G4BestUnit", "G4double");
    (*store)["Time"] = G4AttDef("Time", "Time", "Physics", "G4BestUnit", "G4double");
    (*store)["Pos"]  = G4AttDef("Pos", "Position", "Physics", "G4BestUnit", "G4ThreeVector");
  }
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4AttValue>* Par04Hit::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;
  values->push_back(G4AttValue("HitType", "HadPar04Hit", ""));
  values->push_back(G4AttValue("Energy", G4BestUnit(fEdep, "Energy"), ""));
  values->push_back(G4AttValue("Time", G4BestUnit(fTime, "Time"), ""));
  values->push_back(G4AttValue("Pos", G4BestUnit(fPos, "Length"), ""));
  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04Hit::Print()
{
  std::cout << "\tHit " << fEdep / MeV << " MeV at " << fPos / cm << " cm rotation " << fRot
            << " (R,phi,z)= (" << fRhoId << ", " << fPhiId << ", " << fZId << "), " << fTime << " ns"
            << std::endl;
}
