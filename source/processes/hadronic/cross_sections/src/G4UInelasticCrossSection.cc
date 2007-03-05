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
// 12.08.06 V.Ivanchenko - first implementation
// 22.01.07 V.Ivanchenko - add GetIsoZACrossSection
// 05.03.07 V.Ivanchenko - use G4NucleonNuclearCrossSection
//
//


#include "G4UInelasticCrossSection.hh"

#include "G4ParticleTable.hh"
#include "G4GlauberGribovCrossSection.hh"
#include "G4NucleonNuclearCrossSection.hh"
#include "G4UPiNuclearCrossSection.hh"
#include "G4HadronCrossSections.hh"
#include "G4ParticleDefinition.hh"
#include "G4Element.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UInelasticCrossSection::G4UInelasticCrossSection(const G4ParticleDefinition* p) 
{
  hasGlauber = false;
  thEnergy   = 50.*GeV;
  fGlauber   = new G4GlauberGribovCrossSection();
  fGheisha   = new G4HadronCrossSections();
  fNucleon   = 0;
  fUPi       = 0;
  fGheisha   = 0;
  if(p == G4Proton::Proton() || p == G4Neutron::Neutron())
    fNucleon = new G4NucleonNuclearCrossSection();
  else if(p == G4PionPlus::PionPlus() || p == G4PionMinus::PionMinus())
    fUPi = new G4UPiNuclearCrossSection();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UInelasticCrossSection::~G4UInelasticCrossSection()
{
  delete fGlauber;
  if(fNucleon) delete fNucleon;
  if(fUPi)     delete fUPi;
  if(fGheisha) delete fGheisha;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4UInelasticCrossSection::IsApplicable(const G4DynamicParticle* dp, 
					    const G4Element*  elm)
{
  return IsZAApplicable(dp, elm->GetZ(), elm->GetN());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4UInelasticCrossSection::IsZAApplicable(const G4DynamicParticle* dp, 
					      G4double Z, G4double A)
{
  G4bool res = false;
  if(fNucleon || fUPi) res = true;
  else res = fGheisha->IsApplicable(dp, Z, A);
  if(verboseLevel > 1) 
    G4cout << "G4UInelasticCrossSection::IsApplicable  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()
	   << G4endl;
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UInelasticCrossSection::GetCrossSection(const G4DynamicParticle* dp, 
						   const G4Element* elm, 
						   G4double temp)
{
  return GetIsoZACrossSection(dp, elm->GetZ(), elm->GetN(), temp);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UInelasticCrossSection::GetIsoZACrossSection(const G4DynamicParticle* dp, 
							G4double Z,
							G4double A, 
							G4double)
{
  G4double cross = 0.0;
  G4double ekin = dp->GetKineticEnergy();
  G4int iz = G4int(Z + 0.5);
  if(iz > 92) iz = 92;

    // proton and neutron
  if(fNucleon) { 
    if(iz == 1) cross = fGheisha->GetInelasticCrossSection(dp, Z, A);
    else if(ekin > thEnergy) {
      cross = fGlauber->GetIsoZACrossSection(dp, Z, A);
      cross = theFac[iz]*fGlauber->GetInelasticGlauberGribovXsc();
    } else {
      cross = fNucleon->GetIsoZACrossSection(dp, Z, A);
    }

    // pions
  } else if(fUPi) {
    if(iz == 1) cross = fGheisha->GetInelasticCrossSection(dp, Z, A);
    else if(ekin > thEnergy) {
      cross = fGlauber->GetIsoZACrossSection(dp, Z, A);
      cross = theFac[iz]*fGlauber->GetInelasticGlauberGribovXsc();
    } else {
      cross = fUPi->GetIsoZACrossSection(dp, Z, A);
    }

    //others
  } else {
    if(hasGlauber && ekin > thEnergy) {
      cross = fGlauber->GetIsoZACrossSection(dp, Z, A);
      cross = theFac[iz]*fGlauber->GetInelasticGlauberGribovXsc();
    } else {
      cross = fGheisha->GetInelasticCrossSection(dp, Z, A);
    }
  }

  if(verboseLevel > 1) 
    G4cout << "G4UInelasticCrossSection::GetCrossSection  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()
	   << " in nucleus Z= " << Z << "  A= " << A
	   << G4endl;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UInelasticCrossSection::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  G4DynamicParticle dp;
  G4ParticleDefinition* part = const_cast<G4ParticleDefinition*>(&p);
  dp.SetDefinition(part);
  dp.SetKineticEnergy(thEnergy);
  if(fGlauber->IsZAApplicable(&dp, 2.0, 4.0)) {
    hasGlauber = true;

    G4double A[92] = 
      {
	1.0001, 4.0000, 6.9241, 9.0000, 10.801, 12.011, 14.004, 16.004, 19.000, 20.188,
	23.000, 24.320, 27.000, 28.109, 31.000, 32.094, 35.484, 39.985, 39.135, 40.116,
	45.000, 47.918, 50.998, 52.055, 55.000, 55.910, 59.000, 58.760, 63.617, 65.468,
	69.798, 72.691, 75.000, 79.042, 79.986, 83.887, 85.557, 87.710, 89.000, 91.318,
	93.000, 96.025, 98.000, 101.16, 103.00, 106.51, 107.96, 112.51, 114.91, 118.81,
	121.86, 127.70, 127.00, 131.39, 133.00, 137.42, 139.00, 140.21, 141.00, 144.32,
	145.00, 150.45, 152.04, 157.33, 159.00, 162.57, 165.00, 167.32, 169.00, 173.10,
	175.03, 178.54, 181.00, 183.89, 186.25, 190.27, 192.25, 195.11, 197.00, 200.63,
	204.41, 207.24, 209.00, 209.00, 210.00, 222.00, 223.00, 226.00, 227.00, 232.00,
	231.00, 237.98
      };			 

    G4double csup, csdn;
    for(G4int z=2; z<92; z++) {

      G4double Z = G4double(z);
      csup = fGlauber->GetIsoZACrossSection(&dp, Z, A[z]);
      csup = fGlauber->GetInelasticGlauberGribovXsc();

      // proton and neutron
      if(fNucleon) { 
	csdn = fNucleon->GetIsoZACrossSection(&dp, Z, A[z]);

	// pions
      } else if(fUPi) {
	csdn = fUPi->GetIsoZACrossSection(&dp, Z, A[z]);

	// other
      } else {
	csdn = fGheisha->GetInelasticCrossSection(&dp, Z, A[z]);
      }
      theFac[z] = csdn/csup;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UInelasticCrossSection::DumpPhysicsTable(const G4ParticleDefinition&) 
{
  G4cout << "G4UInelasticCrossSection:"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


