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
// $Id: BrachyPhantomHit.cc,v 1.6 2005/11/22 12:47:35 guatelli Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
#include "BrachyPhantomHit.hh"
#include "G4ios.hh"

G4Allocator<BrachyPhantomHit> BrachyPhantomHitAllocator;

BrachyPhantomHit::BrachyPhantomHit()
{
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
}

const BrachyPhantomHit& BrachyPhantomHit::operator=(const BrachyPhantomHit &right)
{
 xHitPosition = right.xHitPosition;
 zHitPosition = right.zHitPosition;
 yHitPosition = right.yHitPosition;
 energyDeposit = right.energyDeposit;
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
