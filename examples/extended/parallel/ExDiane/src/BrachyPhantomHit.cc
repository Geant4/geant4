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
// Code developed by: S.Guatelli
//
//    ********************************
//    *                              *
//    *     BrachyPhantomHit.cc     *
//    *                              *
//    ********************************
//
// $Id: BrachyPhantomHit.cc,v 1.3 2006/06/29 17:33:23 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
#include "BrachyPhantomHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

G4Allocator<BrachyPhantomHit> BrachyPhantomHitAllocator;

BrachyPhantomHit::BrachyPhantomHit(G4LogicalVolume* logVol,G4int XID,G4int YID,G4int ZID): logicalVolume(logVol),
       xHitPosition(XID),
       zHitPosition(ZID),
       yHitPosition(YID)
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
