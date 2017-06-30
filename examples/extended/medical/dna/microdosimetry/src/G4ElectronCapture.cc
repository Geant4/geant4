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
//---------------------------------------------------------------------------
//
// ClassName:   G4ElectronCapture
//
// Description: The process to kill particles to save CPU
//
// Author:      V.Ivanchenko 31 August 2010
//
//----------------------------------------------------------------------------
//
// $ID$
/// \file G4ElectronCapture.cc
/// \brief Implementation of the G4ElectronCapture class

#include "G4ElectronCapture.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ElectronCapture::G4ElectronCapture(const G4String& regName, G4double ekinlim)
: G4VDiscreteProcess("eCapture", fElectromagnetic),
  fKinEnergyThreshold(ekinlim),
  fRegionName(regName), fpRegion(0)
{
  if(regName == "" || regName == "world") { 
    fRegionName = "DefaultRegionForTheWorld";
  }
  pParticleChange = &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ElectronCapture::~G4ElectronCapture() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ElectronCapture::SetKinEnergyLimit(G4double val)
{
  fKinEnergyThreshold = val;
  if(verboseLevel > 0) {
    G4cout << "### G4ElectronCapture: Tracking cut E(MeV) = " 
        << fKinEnergyThreshold/MeV << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ElectronCapture::BuildPhysicsTable(const G4ParticleDefinition&)
{
  fpRegion = (G4RegionStore::GetInstance())->GetRegion(fRegionName);
  if(fpRegion && verboseLevel > 0) {
    G4cout << "### G4ElectronCapture: Tracking cut E(MeV) = " 
        << fKinEnergyThreshold/MeV << " is assigned to " << fRegionName
        << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4ElectronCapture::IsApplicable(const G4ParticleDefinition&)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4ElectronCapture::PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
    G4double,
    G4ForceCondition* condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;

  G4double limit = DBL_MAX; 
  if(fpRegion) {
    if(aTrack.GetVolume()->GetLogicalVolume()->GetRegion() == fpRegion && 
        aTrack.GetKineticEnergy() < fKinEnergyThreshold) { limit = 0.0; }
  }
  return limit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* G4ElectronCapture::PostStepDoIt(const G4Track& aTrack, 
    const G4Step&)
{
  pParticleChange->Initialize(aTrack);
  pParticleChange->ProposeTrackStatus(fStopAndKill);
  pParticleChange->ProposeLocalEnergyDeposit(aTrack.GetKineticEnergy());
  fParticleChange.SetProposedKineticEnergy(0.0);
  return pParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ElectronCapture::GetMeanFreePath(const G4Track&,G4double,
    G4ForceCondition*)
{
  return DBL_MAX;
}

