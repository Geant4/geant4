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
// Code developed by:
// S. Agostinelli, F. Foppiano, S. Garelli , M. Tropeano, S.Guatelli
//
//    ********************************
//    *                              *
//    *     BrachyPhantomHit.cc     *
//    *                              *
//    ********************************
//
// $Id: BrachyPhantomHit.cc,v 1.5 2003/05/27 08:37:55 guatelli Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
#include "BrachyPhantomHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

G4Allocator<BrachyPhantomHit> BrachyPhantomHitAllocator;

BrachyPhantomHit::BrachyPhantomHit(G4LogicalVolume* logVol,G4int XID,G4int YID,G4int ZID)
  :logicalVolume(logVol),xHitPosition(XID),zHitPosition(ZID),yHitPosition(YID)
{
 energyDeposit = 0;
}

BrachyPhantomHit::~BrachyPhantomHit()
{
}

BrachyPhantomHit::BrachyPhantomHit(const BrachyPhantomHit &right)
  : G4VHit()
{
 xHitPosition = right.xHitPosition;
 zHitPosition = right.zHitPosition;
 yHitPosition = right.yHitPosition;
 energyDeposit = right.energyDeposit;
 hitPosition = right.hitPosition;
 rotation = right.rotation;
 logicalVolume = right.logicalVolume;
}

const BrachyPhantomHit& BrachyPhantomHit::operator=(const BrachyPhantomHit &right)
{
 xHitPosition = right.xHitPosition;
 zHitPosition = right.zHitPosition;
 yHitPosition = right.yHitPosition;
 energyDeposit = right.energyDeposit;
 hitPosition = right.hitPosition;
 rotation = right.rotation;
 logicalVolume = right.logicalVolume;
 return *this;
}

int BrachyPhantomHit::operator==(const BrachyPhantomHit &right) const
{
 return((xHitPosition==right.xHitPosition)&&(zHitPosition==right.zHitPosition)&&(yHitPosition==right.yHitPosition));
}

void BrachyPhantomHit::Draw()
{
}

void BrachyPhantomHit::Print()
{
}
