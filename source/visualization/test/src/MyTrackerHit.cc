// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyTrackerHit.cc,v 1.1 1999-04-16 10:32:37 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MyTrackerHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<MyTrackerHit> MyTrackerHitAllocator;

MyTrackerHit::MyTrackerHit()
{;}

MyTrackerHit::~MyTrackerHit()
{;}

MyTrackerHit::MyTrackerHit(const MyTrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
}

const MyTrackerHit& MyTrackerHit::operator=(const MyTrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
  return *this;
}

int MyTrackerHit::operator==(const MyTrackerHit &right) const
{
  return 0;
}

void MyTrackerHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(0.1);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void MyTrackerHit::Print()
{
}


