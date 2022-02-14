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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file SplitProcess.cc
/// \brief Implementation of the variance reduction particle split class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SplitProcess.hh"
#include "UserTrackInformation.hh"

#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4LogicalVolume.hh"
#include "G4TouchableHandle.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SplitProcess::SplitProcess(G4String regName, G4int nsplit)
:fRegionName(regName), fNSplit(nsplit)
{
    fRegion = G4RegionStore::GetInstance()->FindOrCreateRegion(fRegionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SplitProcess::~SplitProcess() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* SplitProcess::PostStepDoIt(const G4Track& track, const G4Step& step) {
    G4VParticleChange* particleChange(0);
    
    if ( step.GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()
                                                          ->GetRegion() != fRegion ) {
        particleChange = pRegProcess->PostStepDoIt(track, step);
        assert ( 0 != particleChange);
        return particleChange;
    }
    
    if (  fNSplit == 1 ) {
        particleChange = pRegProcess->PostStepDoIt(track, step);
        assert ( 0 != particleChange);
        return particleChange;
    }
    
    UserTrackInformation* parentInformation = 
                          (UserTrackInformation*)(step.GetTrack()->GetUserInformation());
    G4int initialSplitTrackID = parentInformation->GetSplitTrackID();
    
    if ( initialSplitTrackID > 1 ) {
        particleChange = pRegProcess->PostStepDoIt(track, step);
        assert ( 0 != particleChange);
        return particleChange;
    }
    
    G4double weight = track.GetWeight()/fNSplit;
    G4int splitTrackID = 3;
    
    std::vector<G4Track*> secondaries;
    std::vector<G4int> vSplitTrack;
    
    for ( int i = 0; i < fNSplit; i++ ) {
        particleChange = pRegProcess->PostStepDoIt(track, step);
        assert( 0 !=  particleChange);
        particleChange->SetVerboseLevel(0);
        G4Track* newTrack = new G4Track(*(particleChange->GetSecondary(0)));

        secondaries.push_back( newTrack );
        vSplitTrack.push_back( splitTrackID );
        
        splitTrackID++;
    }
    
    parentInformation->SetSplitTrackID(2);
    
    particleChange->SetNumberOfSecondaries(secondaries.size());
    particleChange->SetSecondaryWeightByProcess(true);
    
    std::vector<G4Track*>::iterator iter = secondaries.begin();
    G4int i = 0;
    while( iter != secondaries.end() ) {
        G4Track* newTrack = *iter;
        newTrack->SetWeight(weight);
        
        UserTrackInformation* secondaryInformation = new UserTrackInformation();
        secondaryInformation->SetSplitTrackID(vSplitTrack[i]);
        newTrack->SetUserInformation(secondaryInformation);
        
        particleChange->AddSecondary(newTrack);

        iter++;
        i++;
    }
    
    return particleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
