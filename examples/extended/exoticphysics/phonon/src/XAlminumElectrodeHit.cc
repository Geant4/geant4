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
/// \file exoticphysics/phonon/src/XAlminumElectrodeHit.cc
/// \brief Implementation of the XAlminumElectrodeHit class
//
// $Id$
//

#include "XAlminumElectrodeHit.hh"

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

G4Allocator<XAlminumElectrodeHit> XAlminumElectrodeHitAllocator;

XAlminumElectrodeHit::XAlminumElectrodeHit()
{
  time = 0.;
  fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XAlminumElectrodeHit::~XAlminumElectrodeHit()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XAlminumElectrodeHit::XAlminumElectrodeHit(const XAlminumElectrodeHit &right)
: G4VHit() {
  time = right.time;
  fEdep = right.fEdep;
  fWorldPos = right.fWorldPos;
  fLocalPos = right.fLocalPos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const XAlminumElectrodeHit& XAlminumElectrodeHit::operator=(const XAlminumElectrodeHit &right)
{
  time = right.time;
  fEdep = right.fEdep;
  fWorldPos = right.fWorldPos;
  fLocalPos = right.fLocalPos;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int XAlminumElectrodeHit::operator==(const XAlminumElectrodeHit &/*right*/) const
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XAlminumElectrodeHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fWorldPos);
    circle.SetScreenSize(15);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(0.65,0.65,0.);
    G4VisAttributes attribs(colour);
    attribs.SetStartTime(time);
    attribs.SetEndTime(time+1*ms);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const std::map<G4String,G4AttDef>* XAlminumElectrodeHit::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("XAlminumElectrodeHit",isNew);
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

std::vector<G4AttValue>* XAlminumElectrodeHit::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back(G4AttValue("HitType","XAlminumElectrodeHit",""));

  values->push_back
    (G4AttValue("Time",G4BestUnit(time,"Time"),""));

  values->push_back
    (G4AttValue("EDep",G4BestUnit(fEdep,"Energy"),""));

  values->push_back
    (G4AttValue("Pos",G4BestUnit(fWorldPos,"Length"),""));

  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XAlminumElectrodeHit::Print()
{
  G4cout << "  time " << time/ns << " (nsec) : at " << fLocalPos
         << "  -- fEdep = " << fEdep/eV << " [eV]" << G4endl;
}


