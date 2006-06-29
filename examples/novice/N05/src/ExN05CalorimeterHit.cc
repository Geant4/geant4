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
// $Id: ExN05CalorimeterHit.cc,v 1.6 2006-06-29 17:53:00 gunter Exp $
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
  : G4VHit()
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

G4int ExN05CalorimeterHit::operator==(const ExN05CalorimeterHit &right) const
{
  return (this==&right) ? 1 : 0;
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


