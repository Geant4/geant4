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
// $Id: Par01CalorimeterHit.cc 77659 2013-11-27 08:57:46Z gcosmo $
//

#include "Par01CalorimeterHit.hh"

#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"

G4ThreadLocal G4Allocator<Par01CalorimeterHit>* Par01CalorimeterHitAllocator=0;

Par01CalorimeterHit::Par01CalorimeterHit()
{fLogV=NULL;}

Par01CalorimeterHit::Par01CalorimeterHit(G4LogicalVolume* logVol)
 :fLogV(logVol)
{;}

Par01CalorimeterHit::~Par01CalorimeterHit()
{;}

Par01CalorimeterHit::Par01CalorimeterHit(const Par01CalorimeterHit &right)
  : G4VHit()
{
  fEdep     = right.fEdep;
  fPosition = right.fPosition;
  fRot      = right.fRot;
  fLogV     = right.fLogV;
}

const Par01CalorimeterHit& Par01CalorimeterHit::operator=(const Par01CalorimeterHit &right)
{
  fEdep     = right.fEdep;
  fPosition = right.fPosition;
  fRot      = right.fRot;
  fLogV     = right.fLogV;
  return *this;
}

G4int Par01CalorimeterHit::operator==(const Par01CalorimeterHit &right) const
{
  return (this==&right) ? 1 : 0;
}

void Par01CalorimeterHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Transform3D trans(fRot,fPosition);
    G4VisAttributes attribs;
    const G4VisAttributes* pVA = fLogV->GetVisAttributes();
    if(pVA) attribs = *pVA;
    G4Colour colour(1.,0.,0.);
    attribs.SetColour(colour);
    attribs.SetForceWireframe(false);
    attribs.SetForceSolid(true);
    pVVisManager->Draw(*fLogV,attribs,trans);
  }
}

void Par01CalorimeterHit::Print()
{
}


