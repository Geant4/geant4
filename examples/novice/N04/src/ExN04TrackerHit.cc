
#include "ExN04TrackerHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"


G4Allocator<ExN04TrackerHit> ExN04TrackerHitAllocator;

ExN04TrackerHit::ExN04TrackerHit()
{;}

ExN04TrackerHit::~ExN04TrackerHit()
{;}

ExN04TrackerHit::ExN04TrackerHit(const ExN04TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
}

const ExN04TrackerHit& ExN04TrackerHit::operator=(const ExN04TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
  return *this;
}

int ExN04TrackerHit::operator==(const ExN04TrackerHit &right) const
{
  return 0;
}

void ExN04TrackerHit::Draw()
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

void ExN04TrackerHit::Print()
{;}


