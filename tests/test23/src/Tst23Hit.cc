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
// $Id: Tst23Hit.cc,v 1.1 2001-12-14 14:53:42 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst23Hit.hh"

#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"

G4Allocator<Tst23Hit> Tst23HitAllocator;

Tst23Hit::Tst23Hit()
{;}

Tst23Hit::Tst23Hit(G4int value)
  :id(value)
{;}

Tst23Hit::~Tst23Hit()
{;}

Tst23Hit::Tst23Hit(const Tst23Hit &right)
{
  *this = right;
}

const Tst23Hit& Tst23Hit::operator=(const Tst23Hit &right)
{
  if (this != &right) {
    edep = right.edep;
    id = right.id;
  }
  return *this;
}

int Tst23Hit::operator==(const Tst23Hit &right) const
{
  return 0;
}

void Tst23Hit::Draw()
{
 
}

void Tst23Hit::Print()
{
}


