
#include "ExN04MuonHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"


G4Allocator<ExN04MuonHit> ExN04MuonHitAllocator;

ExN04MuonHit::ExN04MuonHit()
{;}

ExN04MuonHit::~ExN04MuonHit()
{;}

ExN04MuonHit::ExN04MuonHit(const ExN04MuonHit &right)
{
  edep = right.edep;
  pos = right.pos;
}

const ExN04MuonHit& ExN04MuonHit::operator=(const ExN04MuonHit &right)
{
  edep = right.edep;
  pos = right.pos;
  return *this;
}

int ExN04MuonHit::operator==(const ExN04MuonHit &right) const
{
  return 0;
}

void ExN04MuonHit::Draw()
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

void ExN04MuonHit::Print()
{;}


