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
// $Id: PersEx02TrackerHit.cc,v 1.4 2001/07/11 09:58:15 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

#include "PersEx02TrackerHit.hh"

PersEx02TrackerHit::PersEx02TrackerHit()
{;}

PersEx02TrackerHit::~PersEx02TrackerHit()
{;}

PersEx02TrackerHit::PersEx02TrackerHit(const PersEx02TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
}

const PersEx02TrackerHit& PersEx02TrackerHit::operator=(const PersEx02TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
  return *this;
}

int PersEx02TrackerHit::operator==(const PersEx02TrackerHit &right) const
{
  return 0;
}

void PersEx02TrackerHit::Draw()
{
}

void PersEx02TrackerHit::Print()
{
}


