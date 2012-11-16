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
/// \file exoticphysics/phonon/src/XPhononTrackInformation.cc
/// \brief Implementation of the XPhononTrackInformation class
//
// $Id$
//
#include "XPhononTrackInformation.hh"
#include "G4ios.hh"

G4Allocator<XPhononTrackInformation> aTrackInformationAllocator;

XPhononTrackInformation::XPhononTrackInformation()
{
  fK = G4ThreeVector(0., 0., 0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XPhononTrackInformation::XPhononTrackInformation(const G4Track*)
{
  fK = G4ThreeVector(0.,0.,0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XPhononTrackInformation::XPhononTrackInformation(const XPhononTrackInformation* aTrackInfo)
{
  fK = aTrackInfo->GetK(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XPhononTrackInformation::XPhononTrackInformation(G4ThreeVector kNew){
  fK=kNew;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XPhononTrackInformation::~XPhononTrackInformation()
{ ; }

XPhononTrackInformation& XPhononTrackInformation::operator=(const XPhononTrackInformation& aTrackInfo)
{
  fK=aTrackInfo.fK;

  return *this;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XPhononTrackInformation::Print() const
{
  G4cout<<"\nXPhononTrackInformation::Print: Phonon k-vector is "<<fK;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

