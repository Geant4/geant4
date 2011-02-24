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
// $Id: G4BGGNucleonInelasticXS.cc,v 1.13 2011-01-09 02:37:48 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4BGGNucleonInelasticXS
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 13.03.2007
// Modifications:
//
//
// -------------------------------------------------------------------
//

#include "G4BGGNucleonInelasticXS.hh"
#include "G4GlauberGribovCrossSection.hh"
#include "G4NucleonNuclearCrossSection.hh"
#include "G4HadronNucleonXsc.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4NistManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BGGNucleonInelasticXS::G4BGGNucleonInelasticXS(const G4ParticleDefinition* p)
 : G4VCrossSectionDataSet("Barashenkov-Glauber-Gribov")
{
  verboseLevel = 0;
  fGlauberEnergy = 91.*GeV;
  fLowEnergy = 20.*MeV;
  for (G4int i = 0; i < 93; i++) {
    theGlauberFac[i] = 0.0;
    theCoulombFac[i] = 0.0;
  }
  fNucleon = 0;
  fGlauber = 0;
  fHadron  = 0;
  particle = p;
  isProton = false;
  isInitialized = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BGGNucleonInelasticXS::~G4BGGNucleonInelasticXS()
{
  delete fGlauber;
  delete fNucleon;
  delete fHadron;
}


G4double
G4BGGNucleonInelasticXS::GetZandACrossSection(const G4DynamicParticle* dp,
                                              G4int Z, G4int A, G4double)
{
  G4double cross = 0.0;
  G4double ekin = dp->GetKineticEnergy();
//  G4int iz = G4int(Z);
  if(Z > 92) Z = 92;

  if(ekin <= fLowEnergy) {
    cross = theCoulombFac[Z]*CoulombFactor(ekin, A);
  } else if(Z == 1) {
    if( A < 2) {
      fHadron->GetHadronNucleonXscNS(dp, G4Proton::Proton());
      cross = fHadron->GetInelasticHadronNucleonXsc();
    } else {
      fHadron->GetHadronNucleonXscNS(dp, G4Proton::Proton());
      cross = fHadron->GetInelasticHadronNucleonXsc();
      fHadron->GetHadronNucleonXscNS(dp, G4Neutron::Neutron());
      cross += fHadron->GetInelasticHadronNucleonXsc();
    }
  } else if(ekin > fGlauberEnergy) {
    cross = theGlauberFac[Z]*fGlauber->GetInelasticGlauberGribov(dp, Z, A);
  } else {
    cross = fNucleon->GetZandACrossSection(dp, Z, A);
  }

  if(verboseLevel > 1) 
    G4cout << "G4BGGNucleonInelasticXS::GetCrossSection  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()
	   << " in nucleus Z= " << Z << "  A= " << A
	   << " XS(b)= " << cross/barn 
	   << G4endl;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGNucleonInelasticXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(&p == G4Proton::Proton() || &p == G4Neutron::Neutron()) {
    particle = &p;
    Initialise();
  } else {
    G4cout << "### G4BGGNucleonInelasticXS WARNING: is not applicable to " 
	   << p.GetParticleName()
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGNucleonInelasticXS::DumpPhysicsTable(const G4ParticleDefinition&) 
{
  G4cout << "G4BGGNucleonInelasticXS:"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGNucleonInelasticXS::Initialise() 
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
  G4int A = G4lrint(nist->GetAtomicMassAmu(2));

  G4double csup, csdn;

  if(verboseLevel > 0) G4cout << "### G4BGGNucleonInelasticXS::Initialise for "
			      << particle->GetParticleName() << G4endl;

  for(G4int iz=2; iz<93; iz++) {

    G4double Z = G4double(iz);
    A = G4lrint(nist->GetAtomicMassAmu(iz));

    csup = fGlauber->GetInelasticGlauberGribov(&dp, iz, A);
//    csdn = fNucleon->GetIsoZACrossSection(&dp, Z, A);
    csdn = fNucleon->GetZandACrossSection(&dp, iz, A);

    theGlauberFac[iz] = csdn/csup;
    if(verboseLevel > 0) G4cout << "Z= " << Z <<  "  A= " << A 
				<< " factor= " << theGlauberFac[iz] << G4endl; 
  }
  dp.SetKineticEnergy(fLowEnergy);
  fHadron->GetHadronNucleonXscNS(&dp, G4Proton::Proton());
  theCoulombFac[1] = 
    fHadron->GetInelasticHadronNucleonXsc()/CoulombFactor(fLowEnergy,1);
     
  for(G4int iz=2; iz<93; iz++) {

    G4double Z = G4double(iz);
    A = G4lrint(nist->GetAtomicMassAmu(iz));

    theCoulombFac[iz] = 
      fNucleon->GetZandACrossSection(&dp, iz, A)/CoulombFactor(fLowEnergy,A);

    if(verboseLevel > 0) G4cout << "Z= " << Z <<  "  A= " << A 
				<< " factor= " << theCoulombFac[iz] << G4endl; 
  }
}


G4double G4BGGNucleonInelasticXS::CoulombFactor(G4double kinEnergy, G4int A)
{
  G4double res= 0.0;
  if(kinEnergy <= DBL_MIN) return res;
  else if (A < 2) return kinEnergy*kinEnergy;
  
  G4double elog = std::log10(kinEnergy/GeV);
  G4double aa = A;

  // from G4ProtonInelasticCrossSection
  if(isProton) {

    G4double ff1 = 0.70 - 0.002*aa;           // slope of the drop at medium energies.
    G4double ff2 = 1.00 + 1/aa;               // start of the slope.
    G4double ff3 = 0.8  + 18/aa - 0.002*aa;   // stephight
    res = 1.0 + ff3*(1.0 - (1.0/(1+std::exp(-8*ff1*(elog + 1.37*ff2)))));

    ff1 = 1.   - 1./aa - 0.001*aa; // slope of the rise
    ff2 = 1.17 - 2.7/aa-0.0014*aa; // start of the rise
    res /= (1 + std::exp(-8.*ff1*(elog + 2*ff2)));

  } else {

    // from G4NeutronInelasticCrossSection
    G4double p3 = 0.6 + 13./aa - 0.0005*aa;
    G4double p4 = 7.2449 - 0.018242*aa;
    G4double p5 = 1.36 + 1.8/aa + 0.0005*aa;
    G4double p6 = 1. + 200./aa + 0.02*aa;
    G4double p7 = 3.0 - (aa-70.)*(aa-200.)/11000.;

    G4double firstexp  = std::exp(-p4*(elog + p5));
    G4double secondexp = std::exp(-p6*(elog + p7));

    res = (1.+p3*firstexp/(1. + firstexp))/(1. + secondexp);

  }
  return res;  
}

