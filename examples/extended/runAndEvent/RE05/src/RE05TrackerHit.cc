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
<<<<<<< HEAD
// $Id: RE05TrackerHit.cc 69764 2013-05-14 09:59:36Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
/// \file RE05/src/RE05TrackerHit.cc
/// \brief Implementation of the RE05TrackerHit class
//

#include "RE05TrackerHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"
#include "G4AttValue.hh"
#include "G4AttDef.hh"
#include "G4AttCheck.hh"

G4ThreadLocal G4Allocator<RE05TrackerHit>* RE05TrackerHitAllocator=0;

RE05TrackerHit::RE05TrackerHit()
{;}

RE05TrackerHit::~RE05TrackerHit()
{;}

RE05TrackerHit::RE05TrackerHit(const RE05TrackerHit &right)
  : G4VHit()
{
  edep = right.edep;
  pos = right.pos;
}

const RE05TrackerHit& RE05TrackerHit::operator=(const RE05TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
  return *this;
}

<<<<<<< HEAD
G4int RE05TrackerHit::operator==(const RE05TrackerHit &right) const
=======
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool RE05TrackerHit::operator==(const RE05TrackerHit &right) const
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
{
  return (this==&right) ? true : false;
}

std::map<G4String,G4AttDef> RE05TrackerHit::fAttDefs;

void RE05TrackerHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(0.04);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

const std::map<G4String,G4AttDef>* RE05TrackerHit::GetAttDefs() const
{
  // G4AttDefs have to have long life.  Use static member...
  if (fAttDefs.empty()) {
    fAttDefs["HitType"] =
      G4AttDef("HitType","Type of hit","Physics","","G4String");
  }
  return &fAttDefs;
}

std::vector<G4AttValue>* RE05TrackerHit::CreateAttValues() const
{
  // Create expendable G4AttsValues for picking...
  std::vector<G4AttValue>* attValues = new std::vector<G4AttValue>;
  attValues->push_back
    (G4AttValue("HitType","RE05TrackerHit",""));
  //G4cout << "Checking...\n" << G4AttCheck(attValues, GetAttDefs());
  return attValues;
}

void RE05TrackerHit::Print()
{;}


