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
// $Id: ExN05CalorimeterHit.cc,v 1.3 2001-07-11 09:58:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "ExN05CalorimeterHit.hh"

#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"

G4Allocator<ExN05CalorimeterHit> ExN05CalorimeterHitAllocator;

ExN05CalorimeterHit::ExN05CalorimeterHit()
{pLogV=NULL;}

ExN05CalorimeterHit::ExN05CalorimeterHit(G4LogicalVolume* logVol)
:pLogV(logVol)
{;}

ExN05CalorimeterHit::~ExN05CalorimeterHit()
{;}

ExN05CalorimeterHit::ExN05CalorimeterHit(const ExN05CalorimeterHit &right)
{
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
}

const ExN05CalorimeterHit& ExN05CalorimeterHit::operator=(const ExN05CalorimeterHit &right)
{
  edep = right.edep;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
  return *this;
}

int ExN05CalorimeterHit::operator==(const ExN05CalorimeterHit &right) const
{
  return 0;
}

void ExN05CalorimeterHit::Draw()
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

void ExN05CalorimeterHit::Print()
{
}


