//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: PhotInCalorHit.cc,v 1.4 2006/06/29 16:25:05 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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

