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
/// \file SAXSSensitiveDetectorHit.cc
/// \brief Implementation of the SAXSSensitiveDetectorHit class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "SAXSSensitiveDetectorHit.hh"

#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreadLocal G4Allocator<SAXSSensitiveDetectorHit>* hitAllocator = nullptr; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SAXSSensitiveDetectorHit::SAXSSensitiveDetectorHit() :
  G4VHit()
{
  fTime = 0.;
  fPos = G4ThreeVector(0.,0.,0.);
  fMom = G4ThreeVector(0.,0.,0.);
  fEnergy = 0.;
  fType = -1;
  fTrackID = -1;
  fWeight = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SAXSSensitiveDetectorHit::SAXSSensitiveDetectorHit(const SAXSSensitiveDetectorHit &right):
  G4VHit()
{
  fTrackID = right.fTrackID;
  fTrackIDP = right.fTrackIDP;
  fPos = right.fPos;
  fMom = right.fMom;
  fTime = right.fTime;
  fEnergy = right.fEnergy;
  fType = right.fType;
  fWeight = right.fWeight;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const SAXSSensitiveDetectorHit& SAXSSensitiveDetectorHit::operator=
                                       (const SAXSSensitiveDetectorHit &right)
{
  fTrackID = right.fTrackID;
  fTrackIDP = right.fTrackIDP;
  fPos = right.fPos;
  fMom = right.fMom;
  fTime = right.fTime;
  fEnergy = right.fEnergy;
  fType = right.fType;
  fWeight = right.fWeight;
        
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int SAXSSensitiveDetectorHit::operator==(const SAXSSensitiveDetectorHit &) const
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SAXSSensitiveDetectorHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (pVVisManager)
    {
      G4Circle circle(fPos);
      circle.SetScreenSize(2);
      circle.SetFillStyle(G4Circle::filled);
      G4Colour colour(1.,1.,0.);
      G4VisAttributes attribs(colour);
      circle.SetVisAttributes(attribs);
      pVVisManager->Draw(circle);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SAXSSensitiveDetectorHit::Print()
{
  G4cout << " Hit at time " << fTime/ns
    << " (nsec) - pos(x,y,z) " << fPos/mm
    << " (mm) - mom(x,y,z) " << fMom/eV << " eV" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

