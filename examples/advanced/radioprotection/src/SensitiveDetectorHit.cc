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
// Authors: Susanna Guatelli, susanna@uow.edu.au,
// Authors: Jeremy Davis, jad028@uowmail.edu.au
//
// Code based on basic example B02
//
#include "SensitiveDetectorHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

// THIS IS NECESSARY FOR MT MODE
G4ThreadLocal G4Allocator<SensitiveDetectorHit>* SensitiveDetectorHitAllocator=0;

SensitiveDetectorHit::SensitiveDetectorHit()
 : G4VHit(),
   fEdep(0.)
{}

SensitiveDetectorHit::~SensitiveDetectorHit() {}

SensitiveDetectorHit::SensitiveDetectorHit(const SensitiveDetectorHit& right)
  : G4VHit()
{
  fEdep      = right.fEdep;
}

const SensitiveDetectorHit& SensitiveDetectorHit::operator=(const SensitiveDetectorHit& right)
{
  fEdep      = right.fEdep;

  return *this;
}

G4bool SensitiveDetectorHit::operator==(const SensitiveDetectorHit& right) const
{
  return ( this == &right ) ? true : false;
}

void SensitiveDetectorHit::Draw()
{

}

void SensitiveDetectorHit::Print()
{
  G4cout
     << "HIT: Edep: "
     << std::setw(7) << G4BestUnit(fEdep,"Energy")
      << G4endl;
}

