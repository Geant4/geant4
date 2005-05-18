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
#include "HadrontherapyPhantomHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

G4Allocator<HadrontherapyPhantomHit> HadrontherapyPhantomHitAllocator;

HadrontherapyPhantomHit::HadrontherapyPhantomHit()
{
 energyDeposit = 0;
}

HadrontherapyPhantomHit::~HadrontherapyPhantomHit()
{
}

HadrontherapyPhantomHit::HadrontherapyPhantomHit(const HadrontherapyPhantomHit &right)
  : G4VHit()
{
 xHitID = right.xHitID;
 zHitID = right.zHitID;
 yHitID = right.yHitID;
 energyDeposit = right.energyDeposit;
}

const HadrontherapyPhantomHit& HadrontherapyPhantomHit::operator=(const HadrontherapyPhantomHit &right)
{
 xHitID = right.xHitID;
 zHitID = right.zHitID;
 yHitID = right.yHitID;
 energyDeposit = right.energyDeposit;
 return *this;
}

int HadrontherapyPhantomHit::operator==(const HadrontherapyPhantomHit &right) const
{
 return((xHitID==right.xHitID)&&(zHitID==right.zHitID)&&(yHitID==right.yHitID));
}

void HadrontherapyPhantomHit::Draw()
{
}

void HadrontherapyPhantomHit::Print()
{
}
