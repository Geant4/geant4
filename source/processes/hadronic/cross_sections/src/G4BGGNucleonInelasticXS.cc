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
#include "G4NuclearRadii.hh"

#include "G4CrossSectionDataSetRegistry.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4double llog10 = G4Log(10.);

G4double G4BGGNucleonInelasticXS::theGlauberFacP[93] = {0.0};
G4double G4BGGNucleonInelasticXS::theCoulombFacP[93] = {0.0};
G4double G4BGGNucleonInelasticXS::theGlauberFacN[93] = {0.0};
G4double G4BGGNucleonInelasticXS::theCoulombFacN[93] = {0.0};
G4int    G4BGGNucleonInelasticXS::theA[93] = {0};

#ifdef G4MULTITHREADED
G4Mutex G4BGGNucleonInelasticXS::nucleonInelasticXSMutex = G4MUTEX_INITIALIZER;
#endif

G4BGGNucleonInelasticXS::G4BGGNucleonInelasticXS(const G4ParticleDefinition* p)
 : G4VCrossSectionDataSet("BarashenkovGlauberGribov")
{
  verboseLevel = 0;
  fGlauberEnergy = 91.*GeV;
  fLowEnergy = 14.*MeV;

  fNucleon = nullptr;
  fGlauber = nullptr;
  fHadron  = nullptr;

  theProton= G4Proton::Proton();
  isProton = (theProton == p);
  isMaster = false;
  SetForAllAtomsAndEnergies(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BGGNucleonInelasticXS::~G4BGGNucleonInelasticXS()
{
  delete fHadron;
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
  G4double cross = 0.0;
  G4double ekin = dp->GetKineticEnergy();
  G4int Z = std::min(ZZ, 92);
  if(1 == Z) {
    cross = 1.0115*GetIsoCrossSection(dp,1,1);
  } else {
    if(ekin <= fLowEnergy) {
      cross = (isProton) ? theCoulombFacP[Z] : theCoulombFacN[Z];
      cross *= CoulombFactor(ekin, Z);
    } else if(ekin > fGlauberEnergy) {
      cross = (isProton) ? theGlauberFacP[Z] : theGlauberFacN[Z];
      cross *= fGlauber->GetInelasticGlauberGribov(dp, Z, theA[Z]);
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
  fHadron->HadronNucleonXscNS(dp->GetDefinition(), theProton, 
                           dp->GetKineticEnergy());
  G4double cross = A*fHadron->GetInelasticHadronNucleonXsc();

  if(verboseLevel > 1) {
    G4cout << "G4BGGNucleonInelasticXS::GetIsoCrossSection  for "
           << dp->GetDefinition()->GetParticleName()
           << "  Ekin(GeV)= " << dp->GetKineticEnergy()/CLHEP::GeV
           << " in nucleus Z= " << Z << "  A= " << theA[Z]
           << " XS(b)= " << cross/barn 
           << G4endl;
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGNucleonInelasticXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(fNucleon) { return; } 
  if(&p == theProton || &p == G4Neutron::Neutron()) {
    isProton = (theProton == &p);
  } else {
    G4ExceptionDescription ed;
    ed << "This BGG cross section is applicable only to nucleons and not to " 
       << p.GetParticleName() << G4endl; 
    G4Exception("G4BGGNucleonInelasticXS::BuildPhysicsTable", "had001", 
              FatalException, ed);
    return;
  }

  fNucleon = new G4NucleonNuclearCrossSection();
  fGlauber = new G4ComponentGGHadronNucleusXsc();
  fHadron  = new G4HadronNucleonXsc();

  fNucleon->BuildPhysicsTable(p);

  if(0 == theA[0]) { 
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&nucleonInelasticXSMutex);
    if(0 == theA[0]) { 
#endif
      isMaster = true;
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&nucleonInelasticXSMutex);
#endif
  } else {
    return;
  }

  if(isMaster && 0 == theA[0]) {

    theA[0] = theA[1] = 1;
    G4ThreeVector mom(0.0,0.0,1.0);
    G4DynamicParticle dp(theProton, mom, fGlauberEnergy);

    G4NistManager* nist = G4NistManager::Instance();
    G4double csup, csdn;

    if(verboseLevel > 0) {
      G4cout << "### G4BGGNucleonInelasticXS::Initialise for "
            << p.GetParticleName() << G4endl;
    }
    for(G4int iz=2; iz<93; ++iz) {

      G4int A = G4lrint(nist->GetAtomicMassAmu(iz));
      theA[iz] = A;

      csup = fGlauber->GetInelasticGlauberGribov(&dp, iz, A);
      csdn = fNucleon->GetElementCrossSection(&dp, iz);
      theGlauberFacP[iz] = csdn/csup;
    }

    dp.SetDefinition(G4Neutron::Neutron());
    for(G4int iz=2; iz<93; ++iz) {
      csup = fGlauber->GetInelasticGlauberGribov(&dp, iz, theA[iz]);
      csdn = fNucleon->GetElementCrossSection(&dp, iz);
      theGlauberFacN[iz] = csdn/csup;

      if(verboseLevel > 0) {
       G4cout << "Z= " << iz <<  "  A= " << theA[iz] 
              << " GFactorP= " << theGlauberFacP[iz] 
              << " GFactorN= " << theGlauberFacN[iz] << G4endl; 
      }
    }
    theCoulombFacP[1] = theCoulombFacN[1] = 1.0;
    dp.SetDefinition(theProton);
    dp.SetKineticEnergy(fLowEnergy);
    for(G4int iz=2; iz<93; ++iz) {
      theCoulombFacP[iz] = fNucleon->GetElementCrossSection(&dp, iz)
       /CoulombFactor(fLowEnergy, iz);
    }
    dp.SetDefinition(G4Neutron::Neutron());
    for(G4int iz=2; iz<93; ++iz) {
      theCoulombFacN[iz] = fNucleon->GetElementCrossSection(&dp, iz)
       /CoulombFactor(fLowEnergy, iz);

      if(verboseLevel > 0) {
        G4cout << "Z= " << iz <<  "  A= " << theA[iz] 
               << " CFactorP= " << theCoulombFacP[iz] 
               << " CFactorN= " << theCoulombFacN[iz] << G4endl; 
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BGGNucleonInelasticXS::CoulombFactor(G4double kinEnergy, G4int Z)
{
  G4double res = 0.0;

  if(kinEnergy <= 0.0) { return res; }

  G4double elog = G4Log(kinEnergy/GeV)/llog10;
  G4double aa = theA[Z];
 
  if(isProton) { 

    res = G4NuclearRadii::CoulombFactor(Z, theA[Z], theProton, kinEnergy);

    // from G4ProtonInelasticCrossSection
    if(res > 0.0) {
      G4double ff1 = 5.6  - 0.016*aa;    // slope of the drop at medium energies.
      G4double ff2 = 1.37 + 1.37/aa;     // start of the slope.
      G4double ff3 = 0.8  + 18./aa - 0.002*aa;   // stephight
      res *= (1.0 + ff3*(1.0 - (1.0/(1+G4Exp(-ff1*(elog + ff2))))));
      ff1 = 8.   - 8./aa  - 0.008*aa; // slope of the rise
      ff2 = 2.34 - 5.4/aa - 0.0028*aa; // start of the rise
      res /= (1.0 + G4Exp(-ff1*(elog + ff2)));
    }
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
