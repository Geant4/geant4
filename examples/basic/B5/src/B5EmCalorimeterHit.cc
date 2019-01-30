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
/// \file B5EmCalorimeterHit.cc
/// \brief Implementation of the B5EmCalorimeterHit class

#include "B5EmCalorimeterHit.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4Allocator<B5EmCalorimeterHit>* B5EmCalorimeterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5EmCalorimeterHit::B5EmCalorimeterHit()
: G4VHit(), 
  fCellID(-1), fEdep(0.), fPos(0.), fPLogV(nullptr)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5EmCalorimeterHit::B5EmCalorimeterHit(G4int cellID)
: G4VHit(), 
  fCellID(cellID), fEdep(0.), fPos(0.), fPLogV(nullptr)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5EmCalorimeterHit::~B5EmCalorimeterHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5EmCalorimeterHit::B5EmCalorimeterHit(const B5EmCalorimeterHit &right)
: G4VHit(),
  fCellID(right.fCellID),
  fEdep(right.fEdep),
  fPos(right.fPos),
  fRot(right.fRot),
  fPLogV(right.fPLogV)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const B5EmCalorimeterHit& B5EmCalorimeterHit::operator=(const B5EmCalorimeterHit
                                                          &right)
{
  fCellID = right.fCellID;
  fEdep = right.fEdep;
  fPos = right.fPos;
  fRot = right.fRot;
  fPLogV = right.fPLogV;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B5EmCalorimeterHit::operator==(const B5EmCalorimeterHit &right) const
{
  return (fCellID==right.fCellID);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5EmCalorimeterHit::Draw()
{
  auto visManager = G4VVisManager::GetConcreteInstance();
  if ( ! visManager || (fEdep==0.) ) return;

  // Draw a calorimeter cell with a color corresponding to its energy deposit
  G4Transform3D trans(fRot.inverse(),fPos);
  G4VisAttributes attribs;
  auto pVA = fPLogV->GetVisAttributes();
  if (pVA) attribs = *pVA;
  auto rcol = fEdep/(0.7*GeV);
  if (rcol>1.) rcol = 1.;
  if (rcol<0.4) rcol = 0.4;
  G4Colour colour(rcol,0.,0.);
  attribs.SetColour(colour);
  attribs.SetForceSolid(true);
  visManager->Draw(*fPLogV,attribs,trans);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::map<G4String,G4AttDef>* B5EmCalorimeterHit::GetAttDefs() const
{
  G4bool isNew;
  auto store = G4AttDefStore::GetInstance("B5EmCalorimeterHit",isNew);

  if (isNew) {
      (*store)["HitType"] 
        = G4AttDef("HitType","Hit Type","Physics","","G4String");
      
      (*store)["ID"] 
        = G4AttDef("ID","ID","Physics","","G4int");
      
      (*store)["Energy"] 
        = G4AttDef("Energy", "Energy Deposited", "Physics", "G4BestUnit", 
                   "G4double");
      
      (*store)["Pos"] 
        = G4AttDef("Pos", "Position", "Physics","G4BestUnit", 
                   "G4ThreeVector");
      
      (*store)["LVol"] 
        = G4AttDef("LVol","Logical Volume","Physics","","G4String");
  }
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4AttValue>* B5EmCalorimeterHit::CreateAttValues() const
{
  auto values = new std::vector<G4AttValue>;
  
  values
    ->push_back(G4AttValue("HitType","EmCalorimeterHit",""));
  values
    ->push_back(G4AttValue("ID",G4UIcommand::ConvertToString(fCellID),""));
  values
    ->push_back(G4AttValue("Energy",G4BestUnit(fEdep,"Energy"),""));
  values
    ->push_back(G4AttValue("Pos",G4BestUnit(fPos,"Length"),""));
  
  if (fPLogV)
    values->push_back(G4AttValue("LVol",fPLogV->GetName(),""));
  else
    values->push_back(G4AttValue("LVol"," ",""));
  
  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5EmCalorimeterHit::Print()
{
  G4cout << "  Cell[" << fCellID << "] " << fEdep/MeV << " (MeV)" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
