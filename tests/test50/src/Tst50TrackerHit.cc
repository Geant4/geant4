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
// $Id: Tst50TrackerHit.cc,v 1.6 2006-06-29 22:06:26 gunter Exp $
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
: G4VHit()
{
  edep = right.edep;  
}

const Tst50TrackerHit& Tst50TrackerHit::operator=(const Tst50TrackerHit& right)
{
  edep = right.edep;
  return *this;
}

int Tst50TrackerHit::operator==(const Tst50TrackerHit&) const
{
  return 0;
}

void Tst50TrackerHit::Draw()
{

}

void Tst50TrackerHit::Print()
{
 
}


