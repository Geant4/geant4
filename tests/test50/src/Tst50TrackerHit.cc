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
// $Id: Tst50TrackerHit.cc,v 1.4 2003-05-17 18:11:54 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "Tst50TrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<Tst50TrackerHit> Tst50TrackerHitAllocator;

Tst50TrackerHit::Tst50TrackerHit()
{
  edep = 0.;
}

Tst50TrackerHit::~Tst50TrackerHit() 
{}

Tst50TrackerHit::Tst50TrackerHit(const Tst50TrackerHit& right)
{
  edep      = right.edep;  
}

const Tst50TrackerHit& Tst50TrackerHit::operator=(const Tst50TrackerHit& right)
{
  edep      = right.edep;
  return *this;
}

int Tst50TrackerHit::operator==(const Tst50TrackerHit& right) const
{
  return 0;
}

void Tst50TrackerHit::Draw()
{

}

void Tst50TrackerHit::Print()
{
 
}


