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
// $Id: A01DriftChamberHit.cc,v 1.4 2002-12-13 11:34:33 gunter Exp $
// --------------------------------------------------------------
//
#include "A01DriftChamberHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

G4Allocator<A01DriftChamberHit> A01DriftChamberHitAllocator;

A01DriftChamberHit::A01DriftChamberHit()
{
  layerID = -1;
  time = 0.;
}

A01DriftChamberHit::A01DriftChamberHit(G4int z)
{
  layerID = z;
  time = 0.;
}

A01DriftChamberHit::~A01DriftChamberHit()
{;}

A01DriftChamberHit::A01DriftChamberHit(const A01DriftChamberHit &right)
{
  layerID = right.layerID;
  worldPos = right.worldPos;
  localPos = right.localPos;
  time = right.time;
}

const A01DriftChamberHit& A01DriftChamberHit::operator=(const A01DriftChamberHit &right)
{
  layerID = right.layerID;
  worldPos = right.worldPos;
  localPos = right.localPos;
  time = right.time;
  return *this;
}

int A01DriftChamberHit::operator==(const A01DriftChamberHit &right) const
{
  return 0;
}

void A01DriftChamberHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(worldPos);
    circle.SetScreenSize(2);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,1.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void A01DriftChamberHit::Print()
{
  G4cout << "  Layer[" << layerID << "] : time " << time/ns
         << " (nsec) --- local (x,y) " << localPos.x()
         << ", " << localPos.y() << G4endl;
}


