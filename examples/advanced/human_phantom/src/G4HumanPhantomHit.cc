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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli,University of Wollongong, Australia
// Contributions by F. Ambroglini INFN Perugia, Italy
//

#include "G4HumanPhantomHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"


// THIS IS NECESSARY FOR MT MODE
G4ThreadLocal G4Allocator<G4HumanPhantomHit>* G4HumanPhantomHitAllocator=nullptr;

G4HumanPhantomHit::G4HumanPhantomHit() {}

G4HumanPhantomHit::~G4HumanPhantomHit() {}

G4HumanPhantomHit::G4HumanPhantomHit(const G4HumanPhantomHit& right)
  : G4VHit()
{
  bodyPartID = right.bodyPartID;
  edep = right.edep;
}

const G4HumanPhantomHit& G4HumanPhantomHit::operator=(const G4HumanPhantomHit& right)
{
  bodyPartID = right.bodyPartID;
  edep  = right.edep;
  return *this;
}

G4bool G4HumanPhantomHit::operator==(const G4HumanPhantomHit& right) const
{
  return (this==&right) ? true : false;
}

void G4HumanPhantomHit::Draw()
{
}

void G4HumanPhantomHit::Print()
{
  G4cout << "Energy deposit: " << G4BestUnit(edep,"Energy")
	 << "BodyPartID: " << bodyPartID << G4endl;
}

