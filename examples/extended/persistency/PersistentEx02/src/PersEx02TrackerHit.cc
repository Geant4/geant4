// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PersEx02TrackerHit.cc,v 1.3 1999/11/29 18:33:29 morita Exp $
// GEANT4 tag $Name: geant4-03-00 $
//

#include "PersEx02TrackerHit.hh"

PersEx02TrackerHit::PersEx02TrackerHit()
{;}

PersEx02TrackerHit::~PersEx02TrackerHit()
{;}

PersEx02TrackerHit::PersEx02TrackerHit(const PersEx02TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
}

const PersEx02TrackerHit& PersEx02TrackerHit::operator=(const PersEx02TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
  return *this;
}

int PersEx02TrackerHit::operator==(const PersEx02TrackerHit &right) const
{
  return 0;
}

void PersEx02TrackerHit::Draw()
{
}

void PersEx02TrackerHit::Print()
{
}


