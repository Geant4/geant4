
#include "ExE02TrackerHit.hh"

G4Allocator<ExE02TrackerHit> ExE02TrackerHitAllocator;

ExE02TrackerHit::ExE02TrackerHit()
{;}

ExE02TrackerHit::~ExE02TrackerHit()
{;}

ExE02TrackerHit::ExE02TrackerHit(const ExE02TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
}

const ExE02TrackerHit& ExE02TrackerHit::operator=(const ExE02TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
  return *this;
}

int ExE02TrackerHit::operator==(const ExE02TrackerHit &right) const
{
  return 0;
}

void ExE02TrackerHit::Draw()
{
}

void ExE02TrackerHit::Print()
{
}


