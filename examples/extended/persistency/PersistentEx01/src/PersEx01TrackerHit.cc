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


