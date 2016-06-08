
#include "PersEx01TrackerHit.hh"

G4Allocator<PersEx01TrackerHit> PersEx01TrackerHitAllocator;

PersEx01TrackerHit::PersEx01TrackerHit()
{;}

PersEx01TrackerHit::~PersEx01TrackerHit()
{;}

PersEx01TrackerHit::PersEx01TrackerHit(const PersEx01TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
}

const PersEx01TrackerHit& PersEx01TrackerHit::operator=(const PersEx01TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
  return *this;
}

int PersEx01TrackerHit::operator==(const PersEx01TrackerHit &right) const
{
  return 0;
}

void PersEx01TrackerHit::Draw()
{
}

void PersEx01TrackerHit::Print()
{
}


