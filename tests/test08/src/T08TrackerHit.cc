// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08TrackerHit.cc,v 1.1 1999-01-08 16:35:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "T08TrackerHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"


G4Allocator<T08TrackerHit> T08TrackerHitAllocator;

T08TrackerHit::T08TrackerHit()
{;}

T08TrackerHit::~T08TrackerHit()
{;}

T08TrackerHit::T08TrackerHit(const T08TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
}

const T08TrackerHit& T08TrackerHit::operator=(const T08TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
  return *this;
}

int T08TrackerHit::operator==(const T08TrackerHit &right) const
{
  return 0;
}

void T08TrackerHit::Draw()
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

void T08TrackerHit::Print()
{
}


// This is a forward declarations of an instantiated G4Allocator<Type> object.
// It has been added in order to make code portable for the GNU g++ 
// (release 2.7.2) compiler. 
// Whenever a new Type is instantiated via G4Allocator, it has to be forward
// declared to make object code (compiled with GNU g++) link successfully. 
// 
#ifdef GNU_GCC
  template class G4Allocator<T08TrackerHit>;
#endif
