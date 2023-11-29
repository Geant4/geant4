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
//---------------------------------------------------------------------------
//
// ClassName:   G4LowECapture
//
// Author:      V.Ivanchenko 31 August 2010
//
//----------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4LowECapture.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LowECapture::G4LowECapture(G4double ekinlim)
  : G4VDiscreteProcess("Capture", fElectromagnetic), 
    kinEnergyThreshold(ekinlim), nRegions(0), isIon(false)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LowECapture::~G4LowECapture() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LowECapture::SetKinEnergyLimit(G4double val)
{
  kinEnergyThreshold = val;
  if(verboseLevel > 0) {
    G4cout << "### G4LowECapture: Tracking cut E(MeV) = " 
	   << kinEnergyThreshold/MeV << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LowECapture::AddRegion(const G4String& nam)
{
  G4String r = nam;
  if(r == "" || r == "world" || r == "World") r = "DefaultRegionForTheWorld";
  for(G4int i=0; i<nRegions; ++i) {
    if(regionName[i] == r) { return; }
  } 
  regionName.push_back(r);
  ++nRegions;
  if(verboseLevel > 1) {
    G4cout << "### G4LowECapture: new G4Region <" << r << ">" << G4endl; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LowECapture::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  G4RegionStore* store = G4RegionStore::GetInstance();
  for(G4int i=0; i<nRegions; ++i) {  
    const G4Region* r = store->GetRegion(regionName[i]);
    if(r && verboseLevel > 0) {
      G4cout << "### G4LowECapture: new G4Region <" 
	     << regionName[i] << ">  with tracking cut " 
	     << kinEnergyThreshold/keV << " keV" << G4endl; 
    }
    if(r) { region.push_back(r); }
  }
  nRegions = (G4int)region.size();

  // ions reusing G4GenericIon parameters
  if(part.GetParticleType() == "nucleus") {
    G4String pname = part.GetParticleName();
    if(pname != "deuteron" && pname != "triton" &&
       pname != "alpha"    && pname != "He3" &&
       pname != "alpha+"   && pname != "helium" &&
       pname != "hydrogen") { isIon = true; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4LowECapture::IsApplicable(const G4ParticleDefinition&)
{
  return true;
}

G4double G4LowECapture::PostStepGetPhysicalInteractionLength(
	   const G4Track& aTrack, G4double, G4ForceCondition* condition)
{
  *condition = NotForced;  
  G4double limit = DBL_MAX;
  G4double eLimit = kinEnergyThreshold;
  if(isIon) { 
    eLimit *= aTrack.GetDefinition()->GetPDGMass()/CLHEP::proton_mass_c2; 
  }
  if(aTrack.GetKineticEnergy() < eLimit) { 
    for(G4int i=0; i<nRegions; ++i) {  
      if(aTrack.GetVolume()->GetLogicalVolume()->GetRegion() == region[i]) { 
	limit = 0.0; 
	break;
      }
    }
  }
  return limit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* G4LowECapture::PostStepDoIt(const G4Track& aTrack, 
					       const G4Step&)
{
  pParticleChange->Initialize(aTrack);
  pParticleChange->ProposeTrackStatus(fStopAndKill);
  pParticleChange->ProposeLocalEnergyDeposit(aTrack.GetKineticEnergy());
  return pParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4LowECapture::GetMeanFreePath(const G4Track&, G4double,
					G4ForceCondition*)
{
  return DBL_MAX;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


