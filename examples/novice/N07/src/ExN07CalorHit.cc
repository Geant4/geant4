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
//
// $Id: ExN07CalorHit.cc,v 1.2 2003/05/21 22:01:17 asaim Exp $
// GEANT4 tag $Name: geant4-05-02-patch-01 $
//

#include "ExN07CalorHit.hh"

G4Allocator<ExN07CalorHit> ExN07CalorHitAllocator;

ExN07CalorHit::ExN07CalorHit()
:Edep(0.),TrackLength(0.),nStep(0)
{;}

ExN07CalorHit::~ExN07CalorHit()
{;}

ExN07CalorHit::ExN07CalorHit(const ExN07CalorHit& right):G4VHit()
{
  Edep = right.Edep;
  TrackLength = right.TrackLength;
  nStep= right.nStep;
}

const ExN07CalorHit& ExN07CalorHit::operator=(const ExN07CalorHit& right)
{
  Edep = right.Edep;
  TrackLength = right.TrackLength;
  nStep= right.nStep;
  return *this;
}

int ExN07CalorHit::operator==(const ExN07CalorHit& right) const
{
  return (this==&right);
}

void ExN07CalorHit::Draw()
{;}

void ExN07CalorHit::Print()
{;}

