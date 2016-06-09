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
// 06.03.07 V.Ivanchenko - add Initialise function
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
#include "G4NistManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UInelasticCrossSection::G4UInelasticCrossSection(const G4ParticleDefinition* p) 
{
  hasGlauber = false;
  thEnergy   = 100.*GeV;
  fGlauber   = new G4GlauberGribovCrossSection();
  fGheisha   = G4HadronCrossSections::Instance();
  fNucleon   = 0;
  fUPi       = 0;
  if(p == G4Proton::Proton() || p == G4Neutron::Neutron())
    fNucleon = new G4NucleonNuclearCrossSection();
  else if(p == G4PionPlus::PionPlus() || p == G4PionMinus::PionMinus())
    fUPi = new G4UPiNuclearCrossSection();
  verboseLevel = 0;
  Initialise(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UInelasticCrossSection::~G4UInelasticCrossSection()
{
  delete fGlauber;
  if(fNucleon) delete fNucleon;
  if(fUPi)     delete fUPi;
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
      cross = theFac[iz]*fGlauber->GetInelasticGlauberGribov(dp, Z, A);
    } else {
      cross = fNucleon->GetIsoZACrossSection(dp, Z, A);
    }

    // pions
  } else if(fUPi) {
    if(iz == 1) cross = fGheisha->GetInelasticCrossSection(dp, Z, A);
    else if(ekin > thEnergy) {
      cross = theFac[iz]*fGlauber->GetInelasticGlauberGribov(dp, Z, A);
    } else {
      cross = fUPi->GetInelasticCrossSection(dp, Z, A);
    }

    //others
  } else {
    if(hasGlauber && ekin > thEnergy) {
      cross = theFac[iz]*fGlauber->GetInelasticGlauberGribov(dp, Z, A);
    } else if(fGheisha->IsApplicable(dp, Z, A)){
      cross = fGheisha->GetInelasticCrossSection(dp, Z, A);
    }
  }

  if(verboseLevel > 1) 
    G4cout << "G4UInelasticCrossSection::GetCrossSection  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()
	   << " in nucleus Z= " << Z << "  A= " << A
	   << " XS(b)= " << cross/barn 
	   << G4endl;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UInelasticCrossSection::BuildPhysicsTable(const G4ParticleDefinition&)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UInelasticCrossSection::DumpPhysicsTable(const G4ParticleDefinition&) 
{
  G4cout << "G4UInelasticCrossSection:"<<G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UInelasticCrossSection::Initialise(const G4ParticleDefinition* p) 
{
  G4ParticleDefinition* part = const_cast<G4ParticleDefinition*>(p);
  G4ThreeVector mom(0.0,0.0,1.0);
  G4DynamicParticle dp(part, mom, thEnergy);

  G4NistManager* nist = G4NistManager::Instance();
  G4double A = nist->GetAtomicMassAmu(2);

  if(fGlauber->IsZAApplicable(&dp, 2.0, A)) {
    hasGlauber = true;

    if(verboseLevel > 0) G4cout << "### G4UInelasticCrossSection::Initialise for "
				<< p->GetParticleName() << G4endl;
    G4double csup, csdn;
    for(G4int iz=2; iz<93; iz++) {

      G4double Z = G4double(iz);
      A = nist->GetAtomicMassAmu(iz);
      csup = fGlauber->GetInelasticGlauberGribov(&dp, Z, A);

      // proton and neutron
      if(fNucleon) { 
	csdn = fNucleon->GetIsoZACrossSection(&dp, Z, A);

	// pions
      } else if(fUPi) {
	csdn = fUPi->GetInelasticCrossSection(&dp, Z, A);

	// other
      } else if(fGheisha->IsApplicable(&dp, Z, A)) {
	csdn = fGheisha->GetInelasticCrossSection(&dp, Z, A);

      } else {
        csdn = csup;
      }
      theFac[iz] = csdn/csup;
      if(verboseLevel > 0) G4cout << "Z= " << Z <<  "  A= " << A 
				  << " factor= " << theFac[iz] << G4endl; 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


