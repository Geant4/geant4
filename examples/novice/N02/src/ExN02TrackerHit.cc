// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02TrackerHit.cc,v 1.2 1999-12-15 14:49:22 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "ExN02TrackerHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"


G4Allocator<ExN02TrackerHit> ExN02TrackerHitAllocator;

ExN02TrackerHit::ExN02TrackerHit()
{;}

ExN02TrackerHit::~ExN02TrackerHit()
{;}

ExN02TrackerHit::ExN02TrackerHit(const ExN02TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
}

const ExN02TrackerHit& ExN02TrackerHit::operator=(const ExN02TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
  return *this;
}

int ExN02TrackerHit::operator==(const ExN02TrackerHit &right) const
{
  return 0;
}

void ExN02TrackerHit::Draw()
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

void ExN02TrackerHit::Print()
{
}


// This is a forward declarations of an instantiated G4Allocator<Type> object.
// It has been added in order to make code portable for the GNU g++ 
// (release 2.7.2) compiler. 
// Whenever a new Type is instantiated via G4Allocator, it has to be forward
// declared to make object code (compiled with GNU g++) link successfully. 
// 
#ifdef GNU_GCC
  template class G4Allocator<ExN02TrackerHit>;
#endif
