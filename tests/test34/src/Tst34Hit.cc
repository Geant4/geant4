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

#include "Tst34Hit.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"

G4Allocator<Tst34Hit>Tst34HitAllocator;

Tst34Hit::Tst34Hit()
{pLogV=NULL;}

Tst34Hit::Tst34Hit(G4LogicalVolume* logVol)
:pLogV(logVol)
{;}

Tst34Hit::~Tst34Hit()
{;}

Tst34Hit::Tst34Hit(const Tst34Hit &right)
:G4VHit()
//@@@ Tst34Hit:Is it right with the init?
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
// @@@@ return 0;
	if ((pos==right.pos) &&  (edep == right.edep)) return true;
	else return false;
	
}

void Tst34Hit::Draw()
{
/*
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Transform3D trans(rot,pos);
    G4VisAttributes attribs;
    const G4VisAttributes* pVA = pLogV->GetVisAttributes();
    if(pVA) attribs = *pVA;
    G4Colour colour(1.,0.,0.);
    attribs.SetColour(colour);
    attribs.SetForceWireframe(false);
    attribs.SetForceSolid(true);
    pVVisManager->Draw(*pLogV,attribs,trans);
  }
*/	
}

void Tst34Hit::Print()
{
}









