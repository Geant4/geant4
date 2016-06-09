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
// $Id$
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
#include "G4SystemOfUnits.hh"
#include "G4GlauberGribovCrossSection.hh"
#include "G4NucleonNuclearCrossSection.hh"
#include "G4HadronNucleonXsc.hh"
#include "G4ComponentSAIDTotalXS.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4NistManager.hh"

#include "G4CrossSectionDataSetRegistry.hh"

G4BGGNucleonElasticXS::G4BGGNucleonElasticXS(const G4ParticleDefinition* p)
 : G4VCrossSectionDataSet("Barashenkov-Glauber") 
{
  verboseLevel = 0;
  fGlauberEnergy = 91.*GeV;
  fLowEnergy = 14.*MeV;
  fSAIDHighEnergyLimit = 1.3*GeV;
  for (G4int i = 0; i < 93; ++i) {
    theGlauberFac[i] = 0.0;
    theCoulombFac[i] = 0.0;
    theA[i] = 1;
  }
  fNucleon = 0;
  fGlauber = 0;
  fHadron  = 0;
  fSAID    = 0;
  particle = p;
  theProton= G4Proton::Proton();
  isProton = false;
  isInitialized = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BGGNucleonElasticXS::~G4BGGNucleonElasticXS()
{
  delete fHadron;
  delete fSAID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool 
G4BGGNucleonElasticXS::IsElementApplicable(const G4DynamicParticle*, G4int,
					   const G4Material*)
{
  return true;
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
					      G4int ZZ, const G4Material*)
{
  // this method should be called only for Z > 1

  G4double cross = 0.0;
  G4double ekin = dp->GetKineticEnergy();
  G4int Z = ZZ;
  if(1 == Z) {
    cross = 1.0115*GetIsoCrossSection(dp,1,1);
  } else {
    if(Z > 92) { Z = 92; }

    if(ekin <= fLowEnergy) {
      cross = theCoulombFac[Z];

    } else if(ekin > fGlauberEnergy) {
      cross = theGlauberFac[Z]*fGlauber->GetElasticGlauberGribov(dp, Z, theA[Z]);
    } else {
      cross = fNucleon->GetElasticCrossSection(dp, Z);
    }
  }

  if(verboseLevel > 1) {
    G4cout << "G4BGGNucleonElasticXS::GetElementCrossSection  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()/CLHEP::GeV
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

  if(ekin <= fSAIDHighEnergyLimit) {
    cross = fSAID->GetElasticIsotopeCrossSection(particle, ekin, 1, 1);
  } else {
    fHadron->GetHadronNucleonXscPDG(dp, theProton);
    cross = theCoulombFac[1]*fHadron->GetElasticHadronNucleonXsc();
  } 
  //  } else if(ekin <= 20*GeV) {
  //  fHadron->GetHadronNucleonXscNS(dp, G4Proton::Proton());
  //  cross = theGlauberFac[1]*fHadron->GetElasticHadronNucleonXsc();
  cross *= A;

  if(verboseLevel > 1) {
    G4cout << "G4BGGNucleonElasticXS::GetIsoCrossSection  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()/CLHEP::GeV
	   << " in nucleus Z= " << Z << "  A= " << A
	   << " XS(b)= " << cross/barn 
	   << G4endl;
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGNucleonElasticXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(&p == theProton || &p == G4Neutron::Neutron()) {
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

  fNucleon = (G4NucleonNuclearCrossSection*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NucleonNuclearCrossSection::Default_Name());
  fGlauber = (G4GlauberGribovCrossSection*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4GlauberGribovCrossSection::Default_Name());
    
  fHadron  = new G4HadronNucleonXsc();
  fSAID    = new G4ComponentSAIDTotalXS();
  fNucleon->BuildPhysicsTable(*particle);
  fGlauber->BuildPhysicsTable(*particle);
  if(particle == theProton) { 
    isProton = true; 
    fSAIDHighEnergyLimit = 3*GeV;
  }

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
  dp.SetKineticEnergy(fSAIDHighEnergyLimit);
  fHadron->GetHadronNucleonXscPDG(&dp, theProton);
  theCoulombFac[1] = 
    fSAID->GetElasticIsotopeCrossSection(particle,fSAIDHighEnergyLimit,1,1)
    /fHadron->GetElasticHadronNucleonXsc();
  if(verboseLevel > 0) {
    G4cout << "Z=1   A=1" << " CoulombFactor= " << theCoulombFac[1] << G4endl; 
  }

  dp.SetKineticEnergy(fLowEnergy);
  for(G4int iz=2; iz<93; iz++) {
    theCoulombFac[iz] = fNucleon->GetElasticCrossSection(&dp, iz);
    if(verboseLevel > 0) {
      G4cout << "Z= " << iz <<  "  A= " << theA[iz]
	     << " factor= " << theCoulombFac[iz] << G4endl; 
    }
  }
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
