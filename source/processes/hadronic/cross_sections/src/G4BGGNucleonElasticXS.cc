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
// $Id: G4BGGNucleonElasticXS.cc,v 1.12 2011-01-09 02:37:48 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

G4BGGNucleonElasticXS::G4BGGNucleonElasticXS(const G4ParticleDefinition*)
 : G4VCrossSectionDataSet("Barashenkov-Glauber") 
{
  verboseLevel = 0;
  fGlauberEnergy = 91.*GeV;
  fLowEnergy = 14.*MeV;
  for (G4int i = 0; i < 93; ++i) {
    theGlauberFac[i] = 0.0;
    theCoulombFac[i] = 0.0;
    theA[i] = 1;
  }
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

G4bool 
G4BGGNucleonElasticXS::IsElementApplicable(const G4DynamicParticle*, G4int Z,
					   const G4Material*)
{
  return (1 < Z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4BGGNucleonElasticXS::IsIsoApplicable(const G4DynamicParticle*, 
					      G4int Z, G4int A,  
					      const G4Element*,
					      const G4Material*)
{
  return (1 == Z && 2 >= A);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4BGGNucleonElasticXS::GetElementCrossSection(const G4DynamicParticle* dp,
					      G4int zElement, 
					      const G4Material*)
{
  // this method should be called only for Z > 1

  G4double cross = 0.0;
  G4double ekin = dp->GetKineticEnergy();
  G4int Z = zElement; 
  if(Z > 92) { Z = 92; }
  else if(Z < 2) { Z = 2; }

  if(ekin <= fLowEnergy) {
    cross = theCoulombFac[Z];

  } else if(ekin > fGlauberEnergy) {
    cross = theGlauberFac[Z]*fGlauber->GetElasticGlauberGribov(dp, Z, theA[Z]);
  } else {
    cross = fNucleon->GetElasticCrossSection(dp, Z);
  }

  if(verboseLevel > 1) {
    G4cout << "G4BGGNucleonElasticXS::GetElementCrossSection  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()
	   << " in nucleus Z= " << Z << "  A= " << theA[Z]
	   << " XS(b)= " << cross/barn 
	   << G4endl;
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4BGGNucleonElasticXS::GetIsoCrossSection(const G4DynamicParticle* dp, 
					  G4int Z, G4int A, 
					  const G4Isotope*,
					  const G4Element*,
					  const G4Material*)
{
  // this method should be called only for Z = 1

  G4double cross = 0.0;
  G4double ekin = dp->GetKineticEnergy();

  if(ekin <= fLowEnergy) {
    cross = theCoulombFac[1];

  } else if(ekin <= 20*GeV) {
    fHadron->GetHadronNucleonXscNS(dp, G4Proton::Proton());
    cross = theGlauberFac[1]*fHadron->GetElasticHadronNucleonXsc();
  } else {
    fHadron->GetHadronNucleonXscPDG(dp, G4Proton::Proton());
    cross = fHadron->GetElasticHadronNucleonXsc();
  }
  cross *= A;

  if(verboseLevel > 1) {
    G4cout << "G4BGGNucleonElasticXS::GetIsoCrossSection  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()
	   << " in nucleus Z= " << Z << "  A= " << A
	   << " XS(b)= " << cross/barn 
	   << G4endl;
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGNucleonElasticXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(&p == G4Proton::Proton() || &p == G4Neutron::Neutron()) {
    particle = &p;

  } else {
    G4cout << "### G4BGGNucleonElasticXS WARNING: is not applicable to " 
	   << p.GetParticleName()
	   << G4endl;
    throw G4HadronicException(__FILE__, __LINE__,
	  "G4BGGNucleonElasticXS::BuildPhysicsTable is used for wrong particle");
    return;
  }

  if(isInitialized) { return; }
  isInitialized = true;

  fNucleon = new G4NucleonNuclearCrossSection();
  fGlauber = new G4GlauberGribovCrossSection();
  fHadron  = new G4HadronNucleonXsc();
  fNucleon->BuildPhysicsTable(*particle);
  fGlauber->BuildPhysicsTable(*particle);
  if(particle == G4Proton::Proton()) { isProton = true; }

  G4ParticleDefinition* part = const_cast<G4ParticleDefinition*>(particle);
  G4ThreeVector mom(0.0,0.0,1.0);
  G4DynamicParticle dp(part, mom, fGlauberEnergy);

  G4NistManager* nist = G4NistManager::Instance();

  G4double csup, csdn;
  G4int A;

  if(verboseLevel > 0) {
    G4cout << "### G4BGGNucleonElasticXS::Initialise for "
	   << particle->GetParticleName() << G4endl;
  }

  for(G4int iz=2; iz<93; iz++) {

    A = G4lrint(nist->GetAtomicMassAmu(iz));
    theA[iz] = A;

    csup = fGlauber->GetElasticGlauberGribov(&dp, iz, A);
    csdn = fNucleon->GetElasticCrossSection(&dp, iz);

    theGlauberFac[iz] = csdn/csup;
    if(verboseLevel > 0) { 
      G4cout << "Z= " << iz <<  "  A= " << A 
	     << " factor= " << theGlauberFac[iz] << G4endl;
    } 
  }
  dp.SetKineticEnergy(20*GeV);
  fHadron->GetHadronNucleonXscPDG(&dp, G4Proton::Proton());
  theGlauberFac[1] = fHadron->GetElasticHadronNucleonXsc();
  fHadron->GetHadronNucleonXscNS(&dp, G4Proton::Proton());
  theGlauberFac[1] /= fHadron->GetElasticHadronNucleonXsc();

  dp.SetKineticEnergy(fLowEnergy);
  fHadron->GetHadronNucleonXscNS(&dp, G4Proton::Proton());
  theCoulombFac[1] = theGlauberFac[1]*fHadron->GetElasticHadronNucleonXsc();

  for(G4int iz=2; iz<93; iz++) {
    theCoulombFac[iz] = fNucleon->GetElasticCrossSection(&dp, iz);
    if(verboseLevel > 0) {
      G4cout << "Z= " << iz <<  "  A= " << theA[iz]
	     << " factor= " << theCoulombFac[iz] << G4endl; 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BGGNucleonElasticXS::CoulombFactor(G4double kinEnergy, G4int A)
{
  G4double res = 0.0;
  if(kinEnergy <= DBL_MIN) { return res; }
  else if(A < 2) { return kinEnergy*kinEnergy; }

  G4double elog = std::log10(kinEnergy/GeV);
  G4double aa = A;

  // from G4ProtonInelasticCrossSection
  G4double f1 = 8.0  - 8.0/aa - 0.008*aa;
  G4double f2 = 2.34 - 5.4/aa - 0.0028*aa;

  res = 1.0/(1.0 + std::exp(-f1*(elog + f2)));
 
  f1 = 5.6 - 0.016*aa;
  f2 = 1.37 + 1.37/aa;
  res *= ( 1.0 + (0.8 + 18./aa - 0.002*aa)/(1.0 + std::exp(f1*(elog + f2))));
  return res;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGNucleonElasticXS::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "The Barashenkov-Glauber-Gribov cross section handles elastic\n"
          << "scattering of protons and neutrons from nuclei using the\n"
          << "Barashenkov parameterization below 91 GeV and the Glauber-Gribov\n"
          << "parameterization above 91 GeV. n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
