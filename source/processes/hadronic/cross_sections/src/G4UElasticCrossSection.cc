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

#include "G4UElasticCrossSection.hh"

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


G4UElasticCrossSection::G4UElasticCrossSection(const G4ParticleDefinition*)
 : G4VCrossSectionDataSet("G4UElasticCrossSection")
{
  verboseLevel = 0;
  hasGlauber = false;
  thEnergy   = 90.*GeV;
  for (G4int i = 0; i < 93; i++) theFac[i] = 0.0;
  fGlauber   = new G4GlauberGribovCrossSection();
  fGheisha   = G4HadronCrossSections::Instance();
  fNucleon   = 0;
  fUPi       = 0;
}


G4UElasticCrossSection::~G4UElasticCrossSection()
{
  delete fGlauber;
  if(fNucleon) delete fNucleon;
  if(fUPi)     delete fUPi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4UElasticCrossSection::IsApplicable(const G4DynamicParticle* dp, 
					      const G4Element*  elm)
{
  return IsIsoApplicable(dp, G4lrint(elm->GetZ()), G4lrint(elm->GetN()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// G4bool G4UElasticCrossSection::IsZAApplicable(const G4DynamicParticle* dp, 
// 						G4double Z, G4double A)
G4bool G4UElasticCrossSection::IsIsoApplicable(const G4DynamicParticle* dp,
                                               G4int Z, G4int A)
{
  G4bool res = false;
  if(fNucleon || fUPi) res = true;
  else res = fGheisha->IsApplicable(dp, Z, A);
  if(verboseLevel > 1) 
    G4cout << "G4UElasticCrossSection::IsApplicable  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()
	   << G4endl;
  return res;
}


G4double
G4UElasticCrossSection::GetCrossSection(const G4DynamicParticle* dp, 
                                        const G4Element* elm, G4double temp)
{
  G4int Z = G4lrint(elm->GetZ());
  G4int N = G4lrint(elm->GetN());

  return GetZandACrossSection(dp, Z, N, temp);
}


G4double
G4UElasticCrossSection::GetZandACrossSection(const G4DynamicParticle* dp, 
                                             G4int Z, G4int A, G4double)
{
  G4double cross = 0.0;
  G4double ekin = dp->GetKineticEnergy();
//  G4int iz = G4int(Z);
  if(Z > 92) Z = 92;

    // proton and neutron
  if(fNucleon) { 
    if(Z == 1) cross = fGheisha->GetElasticCrossSection(dp, Z, A);
    else if(ekin > thEnergy) {
      cross = theFac[Z]*fGlauber->GetElasticGlauberGribov(dp, Z, A);
    } else {
      cross = fNucleon->GetElasticCrossSection(dp, Z, A);
    }

    // pions
  } else if(fUPi) {
    if(Z == 1) cross = fGheisha->GetElasticCrossSection(dp, Z, A);
    else if(ekin > thEnergy) {
      cross = theFac[Z]*fGlauber->GetElasticGlauberGribov(dp, Z, A);
    } else {
      cross = fUPi->GetElasticCrossSection(dp, Z, A);
    }

    //others
  } else {
    if(hasGlauber && ekin > thEnergy) {
      cross = theFac[Z]*fGlauber->GetElasticGlauberGribov(dp, Z, A);
    } else if(fGheisha->IsApplicable(dp, Z, A)) {
      cross = fGheisha->GetElasticCrossSection(dp, Z, A);
    }
  }

  if(verboseLevel > 1) 
    G4cout << "G4UElasticCrossSection::GetCrossSection  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()
	   << " in nucleus Z= " << Z << "  A= " << A
	   << " XS(b)= " << cross/barn 
	   << G4endl;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UElasticCrossSection::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(&p == G4Proton::Proton() || &p == G4Neutron::Neutron()) {
    fNucleon = new G4NucleonNuclearCrossSection();
  } else if(&p == G4PionPlus::PionPlus() || &p == G4PionMinus::PionMinus()) {
    fUPi = new G4UPiNuclearCrossSection();
  }
  Initialise(&p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UElasticCrossSection::DumpPhysicsTable(const G4ParticleDefinition&) 
{
  G4cout << "G4UElasticCrossSection:"<<G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UElasticCrossSection::Initialise(const G4ParticleDefinition* p) 
{
  G4ParticleDefinition* part = const_cast<G4ParticleDefinition*>(p);
  G4ThreeVector mom(0.0,0.0,1.0);
  G4DynamicParticle dp(part, mom, thEnergy);

  G4NistManager* nist = G4NistManager::Instance();
  G4int A = G4lrint(nist->GetAtomicMassAmu(2));

  if(fGlauber->IsIsoApplicable(&dp, 2, A)) {
    hasGlauber = true;

    if(verboseLevel > 0) G4cout << "### G4UElasticCrossSection::Initialise for "
				<< p->GetParticleName() << G4endl;
    G4double csup, csdn;
    for(G4int iz=2; iz<93; iz++) {

      G4double Z = G4double(iz);
      A = G4lrint(nist->GetAtomicMassAmu(iz));
      csup = fGlauber->GetElasticGlauberGribov(&dp, iz, A);
      // proton and neutron
      if(fNucleon) { 
	csdn = fNucleon->GetElasticCrossSection(&dp, iz, A);

	// pions
      } else if(fUPi) {
	csdn = fUPi->GetElasticCrossSection(&dp, iz, A);

	// other
      } else if(fGheisha->IsApplicable(&dp, iz, A)) {
	csdn = fGheisha->GetElasticCrossSection(&dp, iz, A);

      } else {
        csdn = csup;
      }
      theFac[iz] = csdn/csup;
      if(verboseLevel > 0) G4cout << "Z= " << Z <<  "  A= " << A 
				  << " factor= " << theFac[iz] << G4endl; 
    }
  }
}


