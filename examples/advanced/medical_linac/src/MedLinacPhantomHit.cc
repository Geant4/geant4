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
// $Id: MedLinacPhantomHit.cc,v 1.5 2006/06/29 16:04:31 gunter Exp $
//
//
// Code developed by: M. Piergentili

//
#include "MedLinacPhantomHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

G4Allocator<MedLinacPhantomHit> MedLinacPhantomHitAllocator;

MedLinacPhantomHit::MedLinacPhantomHit(G4LogicalVolume* logVol,G4int XID,G4int YID,G4int ZID)
  :logicalVolume(logVol),xHitPosition(XID),zHitPosition(ZID),yHitPosition(YID)
{
 energyDeposit = 0;

 //G4cout << "Hit generated1" << G4endl;

}

MedLinacPhantomHit::~MedLinacPhantomHit()
{
}

MedLinacPhantomHit::MedLinacPhantomHit(const MedLinacPhantomHit &right)
  : G4VHit()
{
 xHitPosition = right.xHitPosition;
 zHitPosition = right.zHitPosition;
 yHitPosition = right.yHitPosition;
 energyDeposit = right.energyDeposit;
 //G4cout << "=====energyDeposit in ImrtPhantomHit1 e'" << energyDeposit << G4endl;
 hitPosition = right.hitPosition;
 rotation = right.rotation;
 logicalVolume = right.logicalVolume;

 //G4cout << "Hit generated2" << G4endl;

}

const MedLinacPhantomHit& MedLinacPhantomHit::operator=(const MedLinacPhantomHit &right)
{
 xHitPosition = right.xHitPosition;
 zHitPosition = right.zHitPosition;
 yHitPosition = right.yHitPosition;
 energyDeposit = right.energyDeposit;
 hitPosition = right.hitPosition;
 rotation = right.rotation;
 logicalVolume = right.logicalVolume;
 //G4cout << " EnergyDep in PhantomHit2 e' "<<" " << energyDeposit  <<G4endl; 
 return *this;
}

int MedLinacPhantomHit::operator==(const MedLinacPhantomHit &right) const
{
 return((xHitPosition==right.xHitPosition)&&(zHitPosition==right.zHitPosition)&&(yHitPosition==right.yHitPosition));
}

void MedLinacPhantomHit::Draw()
{
}

void MedLinacPhantomHit::Print()
{
}
