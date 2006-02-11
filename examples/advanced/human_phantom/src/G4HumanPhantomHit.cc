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
  trackID    = right.trackID;
  bodypartID = right.bodypartID;
  bodypartName = right.bodypartName;
  edep       = right.edep;
  pos        = right.pos;
}


const G4HumanPhantomHit& G4HumanPhantomHit::operator=(const G4HumanPhantomHit& right)
{
  trackID    = right.trackID;
  bodypartID = right.bodypartID;
  bodypartName = right.bodypartName;
  edep       = right.edep;
  pos        = right.pos;
  return *this;
}


G4int G4HumanPhantomHit::operator==(const G4HumanPhantomHit& right) const
{
  return (this==&right) ? 1 : 0;
}


void G4HumanPhantomHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(2.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void G4HumanPhantomHit::Print()
{
  G4cout << " TrackID: " << trackID 
         << "\n BodyPartID: " << bodypartID 
	 << " \t-> " << bodypartName
	 << " \t\t-> Energy deposit: " << G4BestUnit(edep,"Energy")
    //   << "\n Position: " << G4BestUnit(pos,"Length") 
	 << G4endl;
}

