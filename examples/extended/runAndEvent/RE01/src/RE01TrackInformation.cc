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
/// \file runAndEvent/RE01/src/RE01TrackInformation.cc
/// \brief Implementation of the RE01TrackInformation class
//
//
//

#include "RE01TrackInformation.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"    

G4ThreadLocal G4Allocator<RE01TrackInformation> *
                                   aTrackInformationAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01TrackInformation::RE01TrackInformation()
  : G4VUserTrackInformation()
{
    fOriginalTrackID = 0;
    fParticleDefinition = 0;
    fOriginalPosition = G4ThreeVector(0.,0.,0.);
    fOriginalMomentum = G4ThreeVector(0.,0.,0.);
    fOriginalEnergy = 0.;
    fOriginalTime = 0.;
    fTrackingStatus = 1;
    fSourceTrackID = -1;
    fSourceDefinition = 0;
    fSourcePosition = G4ThreeVector(0.,0.,0.);
    fSourceMomentum = G4ThreeVector(0.,0.,0.);
    fSourceEnergy = 0.;
    fSourceTime = 0.;
    fSuspendedStepID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01TrackInformation::RE01TrackInformation(const G4Track* aTrack)
  : G4VUserTrackInformation()
{
    fOriginalTrackID = aTrack->GetTrackID();
    fParticleDefinition = aTrack->GetDefinition();
    fOriginalPosition = aTrack->GetPosition();
    fOriginalMomentum = aTrack->GetMomentum();
    fOriginalEnergy = aTrack->GetTotalEnergy();
    fOriginalTime = aTrack->GetGlobalTime();
    fTrackingStatus = 1;
    fSourceTrackID = -1;
    fSourceDefinition = 0;
    fSourcePosition = G4ThreeVector(0.,0.,0.);
    fSourceMomentum = G4ThreeVector(0.,0.,0.);
    fSourceEnergy = 0.;
    fSourceTime = 0.;
    fSuspendedStepID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01TrackInformation
::RE01TrackInformation(const RE01TrackInformation* aTrackInfo)
  : G4VUserTrackInformation()
{
    fOriginalTrackID = aTrackInfo->fOriginalTrackID;
    fParticleDefinition = aTrackInfo->fParticleDefinition;
    fOriginalPosition = aTrackInfo->fOriginalPosition;
    fOriginalMomentum = aTrackInfo->fOriginalMomentum;
    fOriginalEnergy = aTrackInfo->fOriginalEnergy;
    fOriginalTime = aTrackInfo->fOriginalTime;
    fTrackingStatus = aTrackInfo->fTrackingStatus;
    fSourceTrackID = aTrackInfo->fSourceTrackID;
    fSourceDefinition = aTrackInfo->fSourceDefinition;
    fSourcePosition = aTrackInfo->fSourcePosition;
    fSourceMomentum = aTrackInfo->fSourceMomentum;
    fSourceEnergy = aTrackInfo->fSourceEnergy;
    fSourceTime = aTrackInfo->fSourceTime;
    fSuspendedStepID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01TrackInformation::~RE01TrackInformation()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01TrackInformation& RE01TrackInformation
::operator =(const RE01TrackInformation& aTrackInfo)
{
    fOriginalTrackID = aTrackInfo.fOriginalTrackID;
    fParticleDefinition = aTrackInfo.fParticleDefinition;
    fOriginalPosition = aTrackInfo.fOriginalPosition;
    fOriginalMomentum = aTrackInfo.fOriginalMomentum;
    fOriginalEnergy = aTrackInfo.fOriginalEnergy;
    fOriginalTime = aTrackInfo.fOriginalTime;
    fTrackingStatus = aTrackInfo.fTrackingStatus;
    fSourceTrackID = aTrackInfo.fSourceTrackID;
    fSourceDefinition = aTrackInfo.fSourceDefinition;
    fSourcePosition = aTrackInfo.fSourcePosition;
    fSourceMomentum = aTrackInfo.fSourceMomentum;
    fSourceEnergy = aTrackInfo.fSourceEnergy;
    fSourceTime = aTrackInfo.fSourceTime;
    fSuspendedStepID = -1;

    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE01TrackInformation::SetSourceTrackInformation(const G4Track* aTrack)
{
    fSourceTrackID = aTrack->GetTrackID();
    fSourceDefinition = aTrack->GetDefinition();
    fSourcePosition = aTrack->GetPosition();
    fSourceMomentum = aTrack->GetMomentum();
    fSourceEnergy = aTrack->GetTotalEnergy();
    fSourceTime = aTrack->GetGlobalTime();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE01TrackInformation::Print() const
{
    G4cout 
     << "Source track ID " << fSourceTrackID << " (" 
     << fSourceDefinition->GetParticleName() << ","
     << fSourceEnergy/GeV << "[GeV]) at " << fSourcePosition << G4endl;
    G4cout
      << "Original primary track ID " << fOriginalTrackID << " (" 
      << fParticleDefinition->GetParticleName() << ","
     << fOriginalEnergy/GeV << "[GeV])" << G4endl;
}

