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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#include "HadrontherapyDetectorHit.hh"


//**********************MT
G4ThreadLocal G4Allocator<HadrontherapyDetectorHit>* HadrontherapyDetectorHitAllocator=0;
//**********************MT

HadrontherapyDetectorHit::HadrontherapyDetectorHit()
: G4VHit()
{
 energyDeposit = 0;
}

HadrontherapyDetectorHit::~HadrontherapyDetectorHit()
{
}

HadrontherapyDetectorHit::HadrontherapyDetectorHit(const HadrontherapyDetectorHit& right)
  : G4VHit()
{
 xHitID = right.xHitID;
 zHitID = right.zHitID;
 yHitID = right.yHitID;
 energyDeposit = right.energyDeposit;
}

const HadrontherapyDetectorHit& HadrontherapyDetectorHit::operator=(const HadrontherapyDetectorHit& right)
{
 xHitID = right.xHitID;
 zHitID = right.zHitID;
 yHitID = right.yHitID;
 energyDeposit = right.energyDeposit;
 return *this;
}

G4int HadrontherapyDetectorHit::operator==(const HadrontherapyDetectorHit& right) const
{
 return((xHitID==right.xHitID)&&(zHitID==right.zHitID)&&(yHitID==right.yHitID));
}
