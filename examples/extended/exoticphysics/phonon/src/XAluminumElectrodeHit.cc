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
/// \file exoticphysics/phonon/src/XAluminumElectrodeHit.cc
/// \brief Implementation of the XAluminumElectrodeHit class
//
//
// 20141008  Allocators must be thread-local, and must be pointers

#include "XAluminumElectrodeHit.hh"

#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"

G4ThreadLocal G4Allocator<XAluminumElectrodeHit>* XAluminumElectrodeHitAllocator = 0;

XAluminumElectrodeHit::XAluminumElectrodeHit()
{
  fTime = 0.;
  fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XAluminumElectrodeHit::~XAluminumElectrodeHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XAluminumElectrodeHit::XAluminumElectrodeHit(const XAluminumElectrodeHit &right)
: G4VHit() {
  fTime = right.fTime;
  fEdep = right.fEdep;
  fWorldPos = right.fWorldPos;
  fLocalPos = right.fLocalPos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const XAluminumElectrodeHit& XAluminumElectrodeHit::operator=(const XAluminumElectrodeHit &right)
{
  fTime = right.fTime;
  fEdep = right.fEdep;
  fWorldPos = right.fWorldPos;
  fLocalPos = right.fLocalPos;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool XAluminumElectrodeHit::operator==(const XAluminumElectrodeHit &/*right*/) const
{
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XAluminumElectrodeHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fWorldPos);
    circle.SetScreenSize(15);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(0.65,0.65,0.);
    G4VisAttributes attribs(colour);
    attribs.SetStartTime(fTime);
    attribs.SetEndTime(fTime+1*ms);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const std::map<G4String,G4AttDef>* XAluminumElectrodeHit::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("XAluminumElectrodeHit",isNew);
  if (isNew) {
    G4String HitType("HitType");
    (*store)[HitType] = G4AttDef(HitType,"Hit Type","Physics","","G4String");

    G4String Time("Time");
    (*store)[Time] = G4AttDef(Time,"Time","Physics","G4BestUnit","G4double");

    G4String EDep("EDep");
    (*store)[EDep] = G4AttDef(Time,"EDep","Physics","G4BestUnit","G4double");

    G4String Pos("Pos");
    (*store)[Pos] = G4AttDef(Pos, "Position",
                      "Physics","G4BestUnit","G4ThreeVector");
  }
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4AttValue>* XAluminumElectrodeHit::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back(G4AttValue("HitType","XAluminumElectrodeHit",""));

  values->push_back
    (G4AttValue("Time",G4BestUnit(fTime,"Time"),""));

  values->push_back
    (G4AttValue("EDep",G4BestUnit(fEdep,"Energy"),""));

  values->push_back
    (G4AttValue("Pos",G4BestUnit(fWorldPos,"Length"),""));

  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XAluminumElectrodeHit::Print()
{
  G4cout << "  time " << fTime/ns << " (nsec) : at " << fLocalPos
         << "  -- fEdep = " << fEdep/eV << " [eV]" << G4endl;
}


