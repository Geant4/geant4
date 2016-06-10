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
// $Id: G4BGGNucleonInelasticXS.cc 93682 2015-10-28 10:09:49Z gcosmo $
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
#include "G4SystemOfUnits.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4NucleonNuclearCrossSection.hh"
#include "G4HadronNucleonXsc.hh"
#include "G4ComponentSAIDTotalXS.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

#include "G4CrossSectionDataSetRegistry.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4double llog10 = G4Log(10.);

G4BGGNucleonInelasticXS::G4BGGNucleonInelasticXS(const G4ParticleDefinition* p)
 : G4VCrossSectionDataSet("Barashenkov-Glauber")
{
  verboseLevel = 0;
  fGlauberEnergy = 91.*GeV;
  fLowEnergy = 14.*MeV;
  fHighEnergy = 5.*GeV;
  fSAIDHighEnergyLimit = 1.3*GeV;
  fLowestXSection = millibarn;
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

G4BGGNucleonInelasticXS::~G4BGGNucleonInelasticXS()
{
  delete fSAID;
  delete fHadron;
  // The cross section registry will delete fNucleon
  delete fGlauber;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4BGGNucleonInelasticXS::IsElementApplicable(const G4DynamicParticle*, 
						    G4int, const G4Material*)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4BGGNucleonInelasticXS::IsIsoApplicable(const G4DynamicParticle*, 
						G4int Z, G4int,  
						const G4Element*,
						const G4Material*)
{
  return (1 == Z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4BGGNucleonInelasticXS::GetElementCrossSection(const G4DynamicParticle* dp,
						G4int ZZ, const G4Material*)
{
  // this method should be called only for Z > 1

  G4double cross = 0.0;
  G4double ekin = dp->GetKineticEnergy();
  G4int Z = ZZ;
  if(1 == Z) {
    cross = 1.0115*GetIsoCrossSection(dp,1,1);
  } else if(2 == Z) {
    if(ekin > fGlauberEnergy) {
      cross = theGlauberFac[Z]*fGlauber->GetInelasticGlauberGribov(dp, Z, theA[Z]);
    } else {
      cross = fNucleon->GetElementCrossSection(dp, Z);
    }

  } else {
    if(Z > 92) { Z = 92; }

    if(ekin <= fLowEnergy) {
      cross = theCoulombFac[Z]*CoulombFactor(ekin, Z);
    } else if(ekin > fGlauberEnergy) {
      cross = theGlauberFac[Z]*fGlauber->GetInelasticGlauberGribov(dp, Z, theA[Z]);
    } else {
      cross = fNucleon->GetElementCrossSection(dp, Z);
    }
  }

  if(verboseLevel > 1) {
    G4cout << "G4BGGNucleonInelasticXS::GetCrossSection  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()/CLHEP::GeV
	   << " in nucleus Z= " << Z << "  A= " << theA[Z]
	   << " XS(b)= " << cross/barn 
	   << G4endl;
  }
  //AR-18Dec2013  if(cross <= fLowestXSection) { cross = 0.0; }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4BGGNucleonInelasticXS::GetIsoCrossSection(const G4DynamicParticle* dp,
					    G4int Z, G4int A,
					    const G4Isotope*,
					    const G4Element*,
					    const G4Material*)
{
  // this method should be called only for Z = 1

  G4double cross = 0.0;
  G4double ekin = dp->GetKineticEnergy();

  if(ekin <= fSAIDHighEnergyLimit) {
    cross = fSAID->GetInelasticIsotopeCrossSection(particle, ekin, 1, 1);
  } else if(ekin < fHighEnergy) {
    fHadron->GetHadronNucleonXscNS(dp, theProton);
    cross = (theCoulombFac[0]/ekin + 1)*fHadron->GetInelasticHadronNucleonXsc();
  } else {
    fHadron->GetHadronNucleonXscPDG(dp, theProton);
    cross = (theCoulombFac[1]/ekin + 1)*fHadron->GetInelasticHadronNucleonXsc();
  } 
  cross *= A; 

  if(verboseLevel > 1) {
    G4cout << "G4BGGNucleonInelasticXS::GetCrossSection  for "
	   << dp->GetDefinition()->GetParticleName()
	   << "  Ekin(GeV)= " << dp->GetKineticEnergy()/CLHEP::GeV
	   << " in nucleus Z= " << Z << "  A= " << A
	   << " XS(b)= " << cross/barn 
	   << G4endl;
  }
  //AR-18Dec2013  if(cross <= fLowestXSection) { cross = 0.0; }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGNucleonInelasticXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(&p == theProton || &p == G4Neutron::Neutron()) {
    particle = &p;
  } else {
    G4cout << "### G4BGGNucleonInelasticXS WARNING: is not applicable to " 
	   << p.GetParticleName()
	   << G4endl;
    throw G4HadronicException(__FILE__, __LINE__,
	  "G4BGGNucleonElasticXS::BuildPhysicsTable is used for wrong particle");
    return;
  }

  if(isInitialized) { return; }
  isInitialized = true;

  fNucleon = (G4NucleonNuclearCrossSection*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NucleonNuclearCrossSection::Default_Name());
  fGlauber = new G4ComponentGGHadronNucleusXsc();
  fHadron  = new G4HadronNucleonXsc();
  fSAID    = new G4ComponentSAIDTotalXS();

  fNucleon->BuildPhysicsTable(*particle);
  fGlauber->BuildPhysicsTable(*particle);

  if(particle == theProton) { 
    isProton = true; 
    fSAIDHighEnergyLimit = 2*GeV;
    fHighEnergy = 2*GeV;
  }

  G4ThreeVector mom(0.0,0.0,1.0);
  G4DynamicParticle dp(particle, mom, fGlauberEnergy);

  G4NistManager* nist = G4NistManager::Instance();
  G4int A;

  G4double csup, csdn;

  if(verboseLevel > 0) {
    G4cout << "### G4BGGNucleonInelasticXS::Initialise for "
	   << particle->GetParticleName() << G4endl;
  }
  for(G4int iz=2; iz<93; iz++) {

    A = G4lrint(nist->GetAtomicMassAmu(iz));
    theA[iz] = A;

    csup = fGlauber->GetInelasticGlauberGribov(&dp, iz, A);
    csdn = fNucleon->GetElementCrossSection(&dp, iz);

    theGlauberFac[iz] = csdn/csup;
    if(verboseLevel > 0) {
      G4cout << "Z= " << iz <<  "  A= " << A 
	     << " GlauberFactor= " << theGlauberFac[iz] << G4endl; 
    }
  }
  //const G4Material* mat = 0;

  dp.SetKineticEnergy(fSAIDHighEnergyLimit);
  fHadron->GetHadronNucleonXscNS(&dp, theProton);
  theCoulombFac[0] = fSAIDHighEnergyLimit*
    (fSAID->GetInelasticIsotopeCrossSection(particle,fSAIDHighEnergyLimit,1,1)
     /fHadron->GetInelasticHadronNucleonXsc() - 1);
  
  //G4cout << "Z=1 E(GeV)= " << fSAIDHighEnergyLimit/GeV
  //	 << "  xsNS(b)= " << fHadron->GetInelasticHadronNucleonXsc()/barn;  
  fHadron->GetHadronNucleonXscPDG(&dp, theProton);
  //G4cout << "  xsPDG(b)= " << fHadron->GetInelasticHadronNucleonXsc()/barn;
  //G4cout << "  xsSAID(b)= " << fSAID->GetInelasticIsotopeCrossSection(particle,fSAIDHighEnergyLimit,1,1)/barn << G4endl;

  dp.SetKineticEnergy(fHighEnergy);
  fHadron->GetHadronNucleonXscPDG(&dp, theProton);
  G4double x = fHadron->GetInelasticHadronNucleonXsc();

  //G4cout << "Z=1 E(GeV)= " << fHighEnergy/GeV
  //	 << "  xsPDG(b)= " << fHadron->GetInelasticHadronNucleonXsc()/barn;  

  fHadron->GetHadronNucleonXscNS(&dp, theProton);
  theCoulombFac[1] = fHighEnergy*((theCoulombFac[0]/fHighEnergy + 1)
				  *fHadron->GetInelasticHadronNucleonXsc()/x - 1);

  fHadron->GetHadronNucleonXscNS(&dp, theProton);
  //G4cout <<"  xsNS(b)= "<<fHadron->GetInelasticHadronNucleonXsc()/barn<<G4endl;

  if(verboseLevel > 0) {
    G4cout << "Z=1   A=1" << " CoulombFactor[0]= " << theCoulombFac[0]
	   << " CoulombFactor[1]= " << theCoulombFac[1] << G4endl; 
  }
  theCoulombFac[2] = 1.0;
     
  dp.SetKineticEnergy(fLowEnergy);
  for(G4int iz=3; iz<93; iz++) {
    theCoulombFac[iz] = 
      fNucleon->GetElementCrossSection(&dp, iz)/CoulombFactor(fLowEnergy, iz);

    if(verboseLevel > 0) {
      G4cout << "Z= " << iz <<  "  A= " << theA[iz] 
	     << " CoulombFactor= " << theCoulombFac[iz] << G4endl; 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BGGNucleonInelasticXS::CoulombFactor(G4double kinEnergy, G4int Z)
{
  G4double res= 0.0;
  if(kinEnergy <= 0.0) { return res; }
  else if (Z <= 1) { return kinEnergy*kinEnergy; }
  
  G4double elog = G4Log(kinEnergy/GeV)/llog10;
  G4double aa = theA[Z];

  // from G4ProtonInelasticCrossSection
  if(isProton) {

    G4double ff1 = 5.6  - 0.016*aa;    // slope of the drop at medium energies.
    G4double ff2 = 1.37 + 1.37/aa;     // start of the slope.
    G4double ff3 = 0.8  + 18./aa - 0.002*aa;   // stephight
    res = 1.0 + ff3*(1.0 - (1.0/(1+G4Exp(-ff1*(elog + ff2)))));

    ff1 = 8.   - 8./aa  - 0.008*aa; // slope of the rise
    ff2 = 2.34 - 5.4/aa - 0.0028*aa; // start of the rise
    res /= (1.0 + G4Exp(-ff1*(elog + ff2)));

  } else {

    // from G4NeutronInelasticCrossSection
    G4double p3 = 0.6 + 13./aa - 0.0005*aa;
    G4double p4 = 7.2449 - 0.018242*aa;
    G4double p5 = 1.36 + 1.8/aa + 0.0005*aa;
    G4double p6 = 1. + 200./aa + 0.02*aa;
    G4double p7 = 3.0 - (aa-70.)*(aa-200.)/11000.;

    G4double firstexp  = G4Exp(-p4*(elog + p5));
    G4double secondexp = G4Exp(-p6*(elog + p7));

    res = (1.+p3*firstexp/(1. + firstexp))/(1. + secondexp);

  }
  return res;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGNucleonInelasticXS::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "The Barashenkov-Glauber-Gribov cross section calculates inelastic\n"
          << "scattering of protons and neutrons from nuclei using the\n"
          << "Barashenkov parameterization below 91 GeV and the Glauber-Gribov\n"
          << "parameterization above 91 GeV.  It uses the G4HadronNucleonXsc\n"
          << "cross section component for hydrogen targets, and the\n"
          << "G4ComponentGGHadronNucleusXsc component for other targets.\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
