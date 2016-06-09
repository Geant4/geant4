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
// $Id: G4HumanPhantomHit.cc,v 1.11 2007-05-15 14:45:35 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4HumanPhantomHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<G4HumanPhantomHit> G4HumanPhantomHitAllocator;


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

G4int G4HumanPhantomHit::operator==(const G4HumanPhantomHit& right) const
{
  return (this==&right) ? 1 : 0;
}

void G4HumanPhantomHit::Draw()
{
}

void G4HumanPhantomHit::Print()
{
    G4cout << "Energy deposit: " << G4BestUnit(edep,"Energy")
	 << "BodyPartID: " << bodyPartID << G4endl;
}

