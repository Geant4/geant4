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
// $Id: RE01TrackInformation.cc,v 1.1 2004/11/26 07:37:42 asaim Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//

#include "RE01TrackInformation.hh"
#include "G4ios.hh"

G4Allocator<RE01TrackInformation> aTrackInformationAllocator;

RE01TrackInformation::RE01TrackInformation()
{
    originalTrackID = 0;
    particleDefinition = 0;
    originalPosition = G4ThreeVector(0.,0.,0.);
    originalMomentum = G4ThreeVector(0.,0.,0.);
    originalEnergy = 0.;
    originalTime = 0.;
    trackingStatus = 1;
    sourceTrackID = -1;
    sourceTrackID = -1;
    sourceDefinition = 0;
    sourcePosition = G4ThreeVector(0.,0.,0.);
    sourceMomentum = G4ThreeVector(0.,0.,0.);
    sourceEnergy = 0.;
    sourceTime = 0.;
}

RE01TrackInformation::RE01TrackInformation(const G4Track* aTrack)
{
    originalTrackID = aTrack->GetTrackID();
    particleDefinition = aTrack->GetDefinition();
    originalPosition = aTrack->GetPosition();
    originalMomentum = aTrack->GetMomentum();
    originalEnergy = aTrack->GetTotalEnergy();
    originalTime = aTrack->GetGlobalTime();
    trackingStatus = 1;
    sourceTrackID = -1;
    sourceDefinition = 0;
    sourcePosition = G4ThreeVector(0.,0.,0.);
    sourceMomentum = G4ThreeVector(0.,0.,0.);
    sourceEnergy = 0.;
    sourceTime = 0.;
}

RE01TrackInformation::RE01TrackInformation(const RE01TrackInformation* aTrackInfo)
{
    originalTrackID = aTrackInfo->originalTrackID;
    particleDefinition = aTrackInfo->particleDefinition;
    originalPosition = aTrackInfo->originalPosition;
    originalMomentum = aTrackInfo->originalMomentum;
    originalEnergy = aTrackInfo->originalEnergy;
    originalTime = aTrackInfo->originalTime;
    trackingStatus = aTrackInfo->trackingStatus;
    sourceTrackID = aTrackInfo->sourceTrackID;
    sourceDefinition = aTrackInfo->sourceDefinition;
    sourcePosition = aTrackInfo->sourcePosition;
    sourceMomentum = aTrackInfo->sourceMomentum;
    sourceEnergy = aTrackInfo->sourceEnergy;
    sourceTime = aTrackInfo->sourceTime;
}

RE01TrackInformation::~RE01TrackInformation()
{ ; }

RE01TrackInformation& RE01TrackInformation::operator =(const RE01TrackInformation& aTrackInfo)
{
    originalTrackID = aTrackInfo.originalTrackID;
    particleDefinition = aTrackInfo.particleDefinition;
    originalPosition = aTrackInfo.originalPosition;
    originalMomentum = aTrackInfo.originalMomentum;
    originalEnergy = aTrackInfo.originalEnergy;
    originalTime = aTrackInfo.originalTime;
    trackingStatus = aTrackInfo.trackingStatus;
    sourceTrackID = aTrackInfo.sourceTrackID;
    sourceDefinition = aTrackInfo.sourceDefinition;
    sourcePosition = aTrackInfo.sourcePosition;
    sourceMomentum = aTrackInfo.sourceMomentum;
    sourceEnergy = aTrackInfo.sourceEnergy;
    sourceTime = aTrackInfo.sourceTime;

    return *this;
}

void RE01TrackInformation::SetSourceTrackInformation(const G4Track* aTrack)
{
    sourceTrackID = aTrack->GetTrackID();
    sourceDefinition = aTrack->GetDefinition();
    sourcePosition = aTrack->GetPosition();
    sourceMomentum = aTrack->GetMomentum();
    sourceEnergy = aTrack->GetTotalEnergy();
    sourceTime = aTrack->GetGlobalTime();
}

void RE01TrackInformation::Print() const
{
    G4cout 
     << "Source track ID " << sourceTrackID << " (" << sourceDefinition->GetParticleName() << ","
     << sourceEnergy/GeV << "[GeV]) at " << sourcePosition << G4endl;
    G4cout
     << "Original primary track ID " << originalTrackID << " (" << particleDefinition->GetParticleName() << ","
     << originalEnergy/GeV << "[GeV])" << G4endl;
}

