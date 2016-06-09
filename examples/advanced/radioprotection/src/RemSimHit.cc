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

#include "RemSimHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
G4Allocator<RemSimHit> RemSimHitAllocator;

RemSimHit::RemSimHit()
{;}

RemSimHit::~RemSimHit()
{;}

RemSimHit::RemSimHit(const RemSimHit &right)
  : G4VHit()
{
  edep = right.edep;
  zID = right.zID;
  position = right.position;
}

const RemSimHit& RemSimHit::operator=(const RemSimHit &right)
{
  edep = right.edep;
  zID = right.zID; 
  position = right.position;
  return *this;
}

G4int RemSimHit::operator==(const RemSimHit &right) const
{
  return (this==&right) ? 1 : 0;
}

void RemSimHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(position);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void RemSimHit::Print()
{ 
 G4cout<< 
         "Hit in Astronaut:" 
         << G4BestUnit(edep, "Energy")
	 << " in position: " 
         << G4BestUnit(position, "Length") <<G4endl;
}


