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
//
/// \file B5/src/HadCalorimeterHit.cc
/// \brief Implementation of the B5::HadCalorimeterHit class

#include "HadCalorimeterHit.hh"
#include "DetectorConstruction.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

namespace B5
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4Allocator<HadCalorimeterHit>* HadCalorimeterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadCalorimeterHit::HadCalorimeterHit(G4int columnID,G4int rowID)
: fColumnID(columnID), fRowID(rowID)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool HadCalorimeterHit::operator==(const HadCalorimeterHit &right) const
{
  return ( fColumnID==right.fColumnID && fRowID==right.fRowID );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadCalorimeterHit::Draw()
{
  auto visManager = G4VVisManager::GetConcreteInstance();
  if (! visManager || (fEdep==0.)) return;

  // Draw a calorimeter cell with depth propotional to the energy deposition
  G4Transform3D trans(fRot.inverse(),fPos);
  G4VisAttributes attribs;
  G4Colour colour(1.,0.,0.);
  attribs.SetColour(colour);
  attribs.SetForceSolid(true);
  G4Box box("dummy",15.*cm,15.*cm,1.*m*fEdep/(0.1*GeV));
  visManager->Draw(box,attribs,trans);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::map<G4String,G4AttDef>* HadCalorimeterHit::GetAttDefs() const
{
  G4bool isNew;
  auto store = G4AttDefStore::GetInstance("HadCalorimeterHit",isNew);

  if (isNew) {
    (*store)["HitType"]
      = G4AttDef("HitType","Hit Type","Physics","","G4String");

    (*store)["Column"]
      = G4AttDef("Column","Column ID","Physics","","G4int");

    (*store)["Row"]
      = G4AttDef("Row","Row ID","Physics","","G4int");

    (*store)["Energy"]
      = G4AttDef("Energy","Energy Deposited","Physics","G4BestUnit",
                 "G4double");

    (*store)["Pos"]
      = G4AttDef("Pos", "Position", "Physics","G4BestUnit",
                 "G4ThreeVector");
  }
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4AttValue>* HadCalorimeterHit::CreateAttValues() const
{
  auto values = new std::vector<G4AttValue>;

  values
    ->push_back(G4AttValue("HitType","HadCalorimeterHit",""));
  values
    ->push_back(G4AttValue("Column",G4UIcommand::ConvertToString(fColumnID),
                           ""));
  values
    ->push_back(G4AttValue("Row",G4UIcommand::ConvertToString(fRowID),""));
  values
    ->push_back(G4AttValue("Energy",G4BestUnit(fEdep,"Energy"),""));
  values
    ->push_back(G4AttValue("Pos",G4BestUnit(fPos,"Length"),""));

  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadCalorimeterHit::Print()
{
  G4cout << "  Cell[" << fRowID << ", " << fColumnID << "] "
    << fEdep/MeV << " (MeV) " << fPos << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
