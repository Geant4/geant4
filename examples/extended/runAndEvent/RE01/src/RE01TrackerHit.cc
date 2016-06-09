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
// $Id: RE01TrackerHit.cc,v 1.1 2004/11/26 07:37:42 asaim Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//


#include "RE01TrackerHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"


G4Allocator<RE01TrackerHit> RE01TrackerHitAllocator;

RE01TrackerHit::RE01TrackerHit()
{;}

RE01TrackerHit::~RE01TrackerHit()
{;}

RE01TrackerHit::RE01TrackerHit(const RE01TrackerHit &right)
  : G4VHit()
{
  edep = right.edep;
  pos = right.pos;
  trackID = right.trackID;
}

const RE01TrackerHit& RE01TrackerHit::operator=(const RE01TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
  trackID = right.trackID;
  return *this;
}

G4int RE01TrackerHit::operator==(const RE01TrackerHit &right) const
{
  return (this==&right) ? 1 : 0;
}

void RE01TrackerHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(0.04);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void RE01TrackerHit::Print()
{
  G4cout << "TrackID " << trackID << "   Position " << pos << "       : " << edep/keV << " [keV]" << G4endl;
}


