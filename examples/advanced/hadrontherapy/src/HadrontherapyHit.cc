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
// $Id: HadrontherapyHit.cc,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#include "HadrontherapyHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

G4Allocator<HadrontherapyHit> HadrontherapyHitAllocator;

// -----------------------------------------------------------------
HadrontherapyHit::HadrontherapyHit(G4int j)
  :slice(j)
{
  m_Edep = 0;
}

// ----------------------------------------------------------------
HadrontherapyHit::~HadrontherapyHit()
{
}

// ----------------------------------------------------------------
HadrontherapyHit::HadrontherapyHit(const HadrontherapyHit &right):G4VHit()
{
  slice =  right.slice;
  m_Edep = right.m_Edep;
}

// -----------------------------------------------------------------
const HadrontherapyHit& HadrontherapyHit::operator=(const HadrontherapyHit &right)
{
  m_Edep = right.m_Edep;
  slice =  right.slice;
  return *this;
}

// -----------------------------------------------------------------
int HadrontherapyHit::operator==(const HadrontherapyHit &right) const
{ 
  return (slice==right.slice);
}

// ----------------------------------------------------------------
void HadrontherapyHit::Draw()
{
}

// ----------------------------------------------------------------
void HadrontherapyHit::Print()
{
}
