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
// $Id: PhotInCalorHit.cc,v 1.3 2005-12-09 16:44:21 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

//#define debug

#include "PhotInCalorHit.hh"

G4Allocator<PhotInCalorHit> PhotInCalorHitAllocator;

PhotInCalorHit::PhotInCalorHit(): G4VHit(), Edep(0.), TrackLength(0.), nSteps(0) {}

PhotInCalorHit::~PhotInCalorHit() {}

PhotInCalorHit::PhotInCalorHit(const PhotInCalorHit& right): G4VHit()
{
  Edep        = right.Edep;
  TrackLength = right.TrackLength;
  nSteps      = right.nSteps;
#ifdef debug
  G4cout<<"PhotInCalorHit:init by hit E="<<Edep<<",L="<<TrackLength<<",S="<<nSteps<<G4endl;
#endif
}

const PhotInCalorHit& PhotInCalorHit::operator=(const PhotInCalorHit& right)
{
  Edep        = right.Edep;
  TrackLength = right.TrackLength;
  nSteps      = right.nSteps;
#ifdef debug
  G4cout<<"PhotInCalorHit::init by eq E="<<Edep<<",L="<<TrackLength<<",S="<<nSteps<<G4endl;
#endif
  return *this;
}

int PhotInCalorHit::operator==(const PhotInCalorHit& right) const { return (this==&right);}

void PhotInCalorHit::Draw() {} // User can draw the Hit

void PhotInCalorHit::Print(){} // User can print the Hit

