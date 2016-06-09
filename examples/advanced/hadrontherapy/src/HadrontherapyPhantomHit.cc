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
// $Id: HadrontherapyPhantomHit.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "HadrontherapyPhantomHit.hh"

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


