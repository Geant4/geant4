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

#include "ExGflashHit.hh"

#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"

G4Allocator<ExGflashHit>ExGflashHitAllocator;

ExGflashHit::ExGflashHit()
{pLogV=NULL;}

ExGflashHit::ExGflashHit(G4LogicalVolume* logVol)
:pLogV(logVol)
{;}

ExGflashHit::~ExGflashHit()
{;}

ExGflashHit::ExGflashHit(const ExGflashHit &right)
:G4VHit()
//@@@ ExGflashHit:Is it right with the init?
{
  edep = right.edep;
  pos = right.pos; 
  start =right.start; 
  rot = right.rot;
  pLogV = right.pLogV;
  crystalnumber = right.crystalnumber;
}

const ExGflashHit & ExGflashHit::operator=(const ExGflashHit &right)
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

int ExGflashHit::operator==(const ExGflashHit &right) const
{
// @@@@ return 0;
	if ((pos==right.pos) &&  (edep == right.edep)) return true;
	else return false;
	
}

void ExGflashHit::Draw()
{
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
}

void ExGflashHit::Print()
{
}









