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
// $Id: MedLinacPhantomHit.cc,v 1.4 2005/07/03 23:27:37 mpiergen Exp $
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
