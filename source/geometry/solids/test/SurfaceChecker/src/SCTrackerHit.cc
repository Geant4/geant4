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
// $Id: SCTrackerHit.cc,v 1.1 2005-05-19 13:07:29 link Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SCTrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<SCTrackerHit> SCTrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SCTrackerHit::SCTrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SCTrackerHit::~SCTrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SCTrackerHit::SCTrackerHit(const SCTrackerHit& right)
  : G4VHit()
{
  trackID   = right.trackID;
  chamberNb = right.chamberNb;
  edep      = right.edep;
  pos       = right.pos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const SCTrackerHit& SCTrackerHit::operator=(const SCTrackerHit& right)
{
  trackID   = right.trackID;
  chamberNb = right.chamberNb;
  edep      = right.edep;
  pos       = right.pos;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int SCTrackerHit::operator==(const SCTrackerHit& right) const
{
  return (this==&right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCTrackerHit::Draw()
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCTrackerHit::Print()
{
  G4cout << "  trackID: " << trackID << "  chamberNb: " << chamberNb
         << "  energy deposit: " << G4BestUnit(edep,"Energy")
	 << "  position: " << G4BestUnit(pos,"Length") << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

