
#include "MyTrackerHit.hh"

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
}

void MyTrackerHit::Print()
{
}


