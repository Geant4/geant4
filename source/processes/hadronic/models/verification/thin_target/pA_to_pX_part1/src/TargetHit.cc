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
// $Id: TargetHit.cc,v 1.1 2003-05-27 13:44:49 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "TargetHit.hh"

G4Allocator<TargetHit> TargetHitAllocator;


TargetHit::TargetHit()
{
   cosTheta = 0.;
   Ekin = 0.;
   Charge = 0.;
}


TargetHit::~TargetHit()
{}


TargetHit::TargetHit(const TargetHit& right)
{
  cosTheta = right.cosTheta;
  Ekin = right.Ekin;
  Charge = right.Charge;
}


const TargetHit& TargetHit::operator=(const TargetHit& right)
{
  cosTheta = right.cosTheta;
  Ekin = right.Ekin;
  Charge = right.Charge;
  return *this;
}


int TargetHit::operator==(const TargetHit& right) const
{
  return 0;
}


void TargetHit::Draw()
{}


void TargetHit::Print()
{}


