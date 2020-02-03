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
/// \file TrackingAction.cc
/// \brief Implementation of the TrackingAction class

#include "TrackingAction.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4Region.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4RegionStore.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction() :
    G4UserTrackingAction(),
    RunInitObserver(), fpTargetRegion(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{
    fpTargetRegion = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
    const G4ParticleDefinition* particleDefinition =
        track->GetParticleDefinition();

    if(particleDefinition == G4Electron::Definition() ||
        particleDefinition == G4Gamma::Definition())
    {
        if(fpTargetRegion == 0)
          // target region is initialized after
          //detector construction instantiation
        {
            G4Exception("TrackingAction","fpTargetRegion == 0",
                        FatalException,"No target has been found");
        }

        const G4ThreeVector& position = track->GetPosition();

        int N =  fpTargetRegion->GetNumberOfRootVolumes();
        std::vector<G4LogicalVolume*>::iterator it_logicalVolumeInRegion =
                fpTargetRegion->GetRootLogicalVolumeIterator();

        bool inside_target = false;

        for(int i = 0; i < N ; i++, it_logicalVolumeInRegion++)
        {
            EInside test_status = (*it_logicalVolumeInRegion)->GetSolid()->
                Inside(position) ;
            if(test_status == kInside)
            {
                inside_target = true;
                break;
            }
            /*
            else if (test_status == kSurface)
            {
            }
            */
        }

        if(inside_target == true)
        {
            fNParticleInTarget[particleDefinition]++;
        }
        else
        {
            fNParticleInWorld[particleDefinition]++;
        }
    }
}

void
TrackingAction::Initialize()
{
  fpTargetRegion = G4RegionStore::GetInstance()->GetRegion("Target");
}
