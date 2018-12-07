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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4BGGPionInelasticXS
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 01.10.2003
// Modifications:
//
// -------------------------------------------------------------------
//

#include "G4BGGPionInelasticXS.hh"
#include "G4SystemOfUnits.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4UPiNuclearCrossSection.hh"
#include "G4HadronNucleonXsc.hh"
#include "G4ComponentSAIDTotalXS.hh"

#include "G4Proton.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4NistManager.hh"
#include "G4Pow.hh"

#include "G4HadronicParameters.hh"

G4BGGPionInelasticXS::G4BGGPionInelasticXS(const G4ParticleDefinition* p) 
 : G4VCrossSectionDataSet("Barashenkov-Glauber-Gribov")
{
  verboseLevel = 0;
  fGlauberEnergy = 91.*GeV;
  fLowEnergy = 20.*MeV;
  fSAIDHighEnergyLimit = 2.6*GeV;
  SetMinKinEnergy(0.0);
  SetMaxKinEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );

  for (G4int i = 0; i < 93; i++) {
    theGlauberFac[i] = 0.0;
    theCoulombFac[i] = 0.0;
    theA[i] = 1;
  }
  fPion = nullptr;
  fGlauber = nullptr;
  fHadron  = nullptr;
  fSAID    = nullptr;

  fG4pow   = G4Pow::GetInstance();

  particle = p;
  theProton= G4Proton::Proton();
  isPiplus = (p == G4PionPlus::PionPlus()) ? true : false;
  isInitialized = false;
  SetForAllAtomsAndEnergies(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BGGPionInelasticXS::~G4BGGPionInelasticXS()
{
  delete fHadron;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool 
G4BGGPionInelasticXS::IsElementApplicable(const G4DynamicParticle*, G4int,
					  const G4Material*)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4BGGPionInelasticXS::IsIsoApplicable(const G4DynamicParticle*, 
					     G4int Z, G4int,  
					     const G4Element*,
					     const G4Material*)
{
  return (1 == Z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4BGGPionInelasticXS::GetElementCrossSection(const G4DynamicParticle* dp,
					     G4int ZZ, const G4Material*)
{
  // this method should be called only for Z > 1

  G4double cross = 0.0;
  G4double ekin = dp->GetKineticEnergy();
  G4int Z = std::min(ZZ, 92);
  if(1 == Z) {
    cross = 1.0115*GetIsoCrossSection(dp,1,1);
  } else {
    if(ekin <= fLowEnergy && !isPiplus) {
      cross = theCoulombFac[Z];
    } else if(ekin <= 2*MeV && isPiplus) {
      cross = theCoulombFac[Z]*CoulombFactor(ekin, Z);
    } else if(ekin > fGlauberEnergy) {
      cross = theGlauberFac[Z]*fGlauber->GetInelasticGlauberGribov(dp, Z, theA[Z]);
    } else {
      cross = fPion->GetInelasticCrossSection(dp, Z, theA[Z]);
    }
  }
  if(verboseLevel > 1) {
    G4cout << "G4BGGPionInelasticXS::GetCrossSection  for "
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
G4BGGPionInelasticXS::GetIsoCrossSection(const G4DynamicParticle* dp, 
					 G4int Z, G4int A, 
					 const G4Isotope*,
					 const G4Element*,
					 const G4Material*)
{
  // this method should be called only for Z = 1

  G4double cross = 0.0;
  if(1 == Z) {
    fHadron->GetHadronNucleonXscNS(dp, theProton);
    cross = fHadron->GetInelasticHadronNucleonXsc();
  }
  cross *= A;

  if(verboseLevel > 1) {
    G4cout << "G4BGGPionInelasticXS::GetCrossSection  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()
	   << " in nucleus Z= " << Z << "  A= " << A
	   << " XS(b)= " << cross/barn 
	   << G4endl;
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGPionInelasticXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(&p == G4PionPlus::PionPlus() || &p == G4PionMinus::PionMinus()) {
    particle = &p;
  } else {
    G4cout << "### G4BGGPionInelasticXS WARNING: is not applicable to " 
	   << p.GetParticleName()
	   << G4endl;
    throw G4HadronicException(__FILE__, __LINE__,
	  "G4BGGPionInelasticXS::BuildPhysicsTable is used for wrong particle");
    return;
  }

  if(isInitialized) { return; }
  isInitialized = true;

  fPion    = new G4UPiNuclearCrossSection();
  fGlauber = new G4ComponentGGHadronNucleusXsc();
  fHadron  = new G4HadronNucleonXsc();

  fPion->BuildPhysicsTable(*particle);
  fGlauber->BuildPhysicsTable(*particle);

  G4ThreeVector mom(0.0,0.0,1.0);
  G4DynamicParticle dp(particle, mom, fGlauberEnergy);

  G4NistManager* nist = G4NistManager::Instance();

  G4double csup, csdn;

  if(verboseLevel > 0) {
    G4cout << "### G4BGGPionInelasticXS::Initialise for "
	   << particle->GetParticleName()
	   << " isPiplus: " << isPiplus
	   << G4endl;
  }

  for(G4int iz=2; iz<93; ++iz) {
    G4int A = G4lrint(nist->GetAtomicMassAmu(iz));
    theA[iz] = A;

    csup = fGlauber->GetInelasticGlauberGribov(&dp, iz, A);
    csdn = fPion->GetInelasticCrossSection(&dp, iz, A);

    theGlauberFac[iz] = csdn/csup;
    if(verboseLevel > 0) {
      G4cout << "Z= " << iz <<  "  A= " << A 
	     << " factor= " << theGlauberFac[iz] << G4endl; 
    }
  }
  theCoulombFac[1] = 1.0;
  G4double elim = (isPiplus) ? 2*MeV : fLowEnergy; 
  dp.SetKineticEnergy(elim);

  for(G4int iz=2; iz<93; ++iz) {
    if(isPiplus) { 
      theCoulombFac[iz] = fPion->GetInelasticCrossSection(&dp, iz, theA[iz])
	/CoulombFactor(elim, iz); 
    } else {
      theCoulombFac[iz] = fPion->GetInelasticCrossSection(&dp, iz, theA[iz]);
    }
    if(verboseLevel > 0) {
      G4cout << "Z= " << iz <<  "  A= " << theA[iz] 
	     << " CoulombFactor= " << theCoulombFac[iz] << G4endl; 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BGGPionInelasticXS::CoulombFactor(G4double kinEnergy, G4int Z)
{
  G4int A = theA[Z];
  G4double res= 0.0;
  if(kinEnergy <= DBL_MIN) { return res; }
  else if(A < 2) { return kinEnergy*kinEnergy; }
  
  G4double elog = fG4pow->log10A(6.7*kinEnergy/GeV);
  G4double aa = A;

  // from G4ProtonInelasticCrossSection
  G4double ff1 = 0.70 - 0.002*aa;           // slope of the drop at medium energies.
  G4double ff2 = 1.00 + 1/aa;               // start of the slope.
  G4double ff3 = 0.8  + 18/aa - 0.002*aa;   // stephight
  res = 1.0 + ff3*(1.0 - (1.0/(1+fG4pow->expA(-8*ff1*(elog + 1.37*ff2)))));

  ff1 = 1.   - 1./aa - 0.001*aa; // slope of the rise
  ff2 = 1.17 - 2.7/aa-0.0014*aa; // start of the rise
  res /= (1 + fG4pow->expA(-8.*ff1*(elog + 2*ff2)));
  return res;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
G4BGGPionInelasticXS::CrossSectionDescription(std::ostream& outFile) const 
{
  outFile << "The Barashenkov-Glauber-Gribov cross section handles inelastic\n"
          << "pion scattering from nuclei at all energies.  The Barashenkov\n"
          << "parameterization is used below 91 GeV and the Glauber-Gribov\n"
          << "parameterization is used above 91 GeV.\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
