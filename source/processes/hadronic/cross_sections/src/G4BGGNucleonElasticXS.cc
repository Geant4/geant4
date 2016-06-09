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
// $Id: G4BGGNucleonElasticXS.cc,v 1.3 2008/12/01 16:50:23 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4BGGNucleonElasticXS
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 13.03.2007
// Modifications:
//
//
// -------------------------------------------------------------------
//

#include "G4BGGNucleonElasticXS.hh"
#include "G4GlauberGribovCrossSection.hh"
#include "G4NucleonNuclearCrossSection.hh"
#include "G4HadronNucleonXsc.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4NistManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BGGNucleonElasticXS::G4BGGNucleonElasticXS(const G4ParticleDefinition*) 
{
  verboseLevel = 0;
  fGlauberEnergy = 91.*GeV;
  fLowEnergy = 20.*MeV;
  fNucleon = 0;
  fGlauber = 0;
  fHadron  = 0;
  particle = 0;
  isProton = false;
  isInitialized = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BGGNucleonElasticXS::~G4BGGNucleonElasticXS()
{
  delete fGlauber;
  delete fNucleon;
  delete fHadron;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BGGNucleonElasticXS::GetIsoZACrossSection(const G4DynamicParticle* dp, 
						       G4double Z,
						       G4double A, 
						       G4double)
{
  G4double cross = 0.0;
  G4double ekin = dp->GetKineticEnergy();
  G4int iz = G4int(Z);
  if(iz > 92) iz = 92;

  if(ekin <= fLowEnergy) {
    cross = theCoulombFac[iz];
    if(isProton) { cross *= CoulombFactor(ekin, A); }

  } else if(iz == 1) {
    if( A < 1.5) {
      //fHadron->GetHadronNucleonXscPDG(dp, G4Proton::Proton());
      //fHadron->GetHadronNucleonXscEL(dp, G4Proton::Proton());
      //fHadron->GetHadronNucleonXscNS(dp, G4Proton::Proton());
      fHadron->GetHadronNucleonXscMK(dp, G4Proton::Proton());
      cross = fHadron->GetElasticHadronNucleonXsc();
    } else {
      fHadron->GetHadronNucleonXscMK(dp, G4Proton::Proton());
      cross = fHadron->GetElasticHadronNucleonXsc();
      fHadron->GetHadronNucleonXscMK(dp, G4Neutron::Neutron());
      cross += fHadron->GetElasticHadronNucleonXsc();
    }
  } else if(ekin > fGlauberEnergy) {
    cross = theGlauberFac[iz]*fGlauber->GetElasticGlauberGribov(dp, Z, A);
  } else {
    cross = fNucleon->GetElasticCrossSection(dp, Z, A);
  }

  if(verboseLevel > 1) 
    G4cout << "G4BGGNucleonElasticXS::GetCrossSection  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()
	   << " in nucleus Z= " << Z << "  A= " << A
	   << " XS(b)= " << cross/barn 
	   << G4endl;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGNucleonElasticXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(&p == G4Proton::Proton() || &p == G4Neutron::Neutron()) {
    particle = &p;
    Initialise();
  } else {
    G4cout << "### G4BGGNucleonElasticXS WARNING: is not applicable to " 
	   << p.GetParticleName()
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGNucleonElasticXS::DumpPhysicsTable(const G4ParticleDefinition&) 
{
  G4cout << "G4BGGNucleonElasticXS:"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGNucleonElasticXS::Initialise() 
{
  if(isInitialized) return;
  isInitialized = true;

  fNucleon = new G4NucleonNuclearCrossSection();
  fGlauber = new G4GlauberGribovCrossSection();
  fHadron  = new G4HadronNucleonXsc();
  fNucleon->BuildPhysicsTable(*particle);
  fGlauber->BuildPhysicsTable(*particle);
  if(particle == G4Proton::Proton()) isProton = true;

  G4ParticleDefinition* part = const_cast<G4ParticleDefinition*>(particle);
  G4ThreeVector mom(0.0,0.0,1.0);
  G4DynamicParticle dp(part, mom, fGlauberEnergy);

  G4NistManager* nist = G4NistManager::Instance();

  G4double A, csup, csdn;

  if(verboseLevel > 0) G4cout << "### G4BGGNucleonElasticXS::Initialise for "
			      << particle->GetParticleName() << G4endl;

  for(G4int iz=2; iz<93; iz++) {

    G4double Z = G4double(iz);
    A = nist->GetAtomicMassAmu(iz);

    csup = fGlauber->GetElasticGlauberGribov(&dp, Z, A);
    csdn = fNucleon->GetElasticCrossSection(&dp, Z, A);

    theGlauberFac[iz] = csdn/csup;
    if(verboseLevel > 0) G4cout << "Z= " << Z <<  "  A= " << A 
				<< " factor= " << theGlauberFac[iz] << G4endl; 
  }
  dp.SetKineticEnergy(fLowEnergy);
  fHadron->GetHadronNucleonXscMK(&dp, G4Proton::Proton());
  theCoulombFac[1] = fHadron->GetElasticHadronNucleonXsc();
  if(isProton) { theCoulombFac[1] /= CoulombFactor(fLowEnergy, 1.0); }

  for(G4int iz=2; iz<93; iz++) {

    G4double Z = G4double(iz);
    A = nist->GetAtomicMassAmu(iz);

    theCoulombFac[iz] = fNucleon->GetElasticCrossSection(&dp, Z, A);
    if(isProton) { theCoulombFac[iz] /= CoulombFactor(fLowEnergy, A); }
    if(verboseLevel > 0) G4cout << "Z= " << Z <<  "  A= " << A 
				<< " factor= " << theCoulombFac[iz] << G4endl; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BGGNucleonElasticXS::CoulombFactor(G4double kinEnergy, G4double A)
{
  G4double res= 0.0;
  if(kinEnergy <= DBL_MIN) return res;
  else if(A < 1.5) return kinEnergy*kinEnergy;

  G4double elog = std::log10(kinEnergy/GeV);

  // from G4ProtonInelasticCrossSection
  G4double f1 = 8.0  - 8.0/A - 0.008*A;
  G4double f2 = 2.34 - 5.4/A - 0.0028*A;

  res = 1.0/(1.0 + std::exp(-f1*(elog + f2)));
 
  f1 = 5.6 - 0.016*A;
  f2 = 1.37 + 1.37/A;
  res *= ( 1.0 + (0.8 + 18./A - 0.002*A)/(1.0 + std::exp(f1*(elog + f2))));
  return res;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


