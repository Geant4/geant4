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
// $Id: Tst20CalorHit.cc,v 1.5 2007-11-09 18:33:00 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $


#include "Tst20CalorHit.hh"

G4Allocator<Tst20CalorHit> Tst20CalorHitAllocator;


Tst20CalorHit::Tst20CalorHit()
{
   energyDeposit = 0.; 
   trackLength = 0.;
}


Tst20CalorHit::~Tst20CalorHit()
{ }



Tst20CalorHit::Tst20CalorHit(const Tst20CalorHit& right) : G4VHit()
{
  energyDeposit = right.energyDeposit; 
  trackLength = right.trackLength;
}


const Tst20CalorHit& Tst20CalorHit::operator=(const Tst20CalorHit& right)
{
  energyDeposit = right.energyDeposit; 
  trackLength = right.trackLength;
  return *this;
}

bool Tst20CalorHit::operator==(const Tst20CalorHit& right) const
{
  return (this==&right) ? true : false;
}


void Tst20CalorHit::Print()
{ }


void Tst20CalorHit::AddEnergyDeposit(G4double energy, G4double length) 
{
  energyDeposit += energy; 
  trackLength += length;
}


