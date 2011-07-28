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

#include "Tst34Hit.hh"

#include "G4LogicalVolume.hh"
#include "G4ios.hh"

G4Allocator<Tst34Hit>Tst34HitAllocator;

Tst34Hit::Tst34Hit() : pLogV(0)
{;}

Tst34Hit::Tst34Hit(G4LogicalVolume* logVol) : pLogV(logVol)
{;}

Tst34Hit::~Tst34Hit()
{;}

Tst34Hit::Tst34Hit(const Tst34Hit &right) : G4VHit(right)
{
  edep = right.edep;
  pos = right.pos; 
  start =right.start; 
  rot = right.rot;
  pLogV = right.pLogV;
  crystalnumber = right.crystalnumber;
}

const Tst34Hit & Tst34Hit::operator=(const Tst34Hit &right)
{
  edep = right.edep;
  start =right.start; 
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
  crystalnumber = right.crystalnumber;
  crystalnumber = right.crystalnumber;
  return *this;
}

int Tst34Hit::operator==(const Tst34Hit &right) const
{
  if ((pos==right.pos) && (edep == right.edep)) return true;
  else return false;
  
}

void Tst34Hit::Draw()
{
}

void Tst34Hit::Print()
{
}
