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
// File name:     G4BGGPionElasticXS
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 01.10.2003
// Modifications:
//
// -------------------------------------------------------------------
//

#include "G4BGGPionElasticXS.hh"
#include "G4SystemOfUnits.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4UPiNuclearCrossSection.hh"
#include "G4HadronNucleonXsc.hh"
#include "G4NuclearRadii.hh"

#include "G4Proton.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4NistManager.hh"
#include "G4HadronicParameters.hh"
#include "G4Pow.hh"

G4double G4BGGPionElasticXS::theGlauberFacPiPlus[93] = {0.0};
G4double G4BGGPionElasticXS::theCoulombFacPiPlus[93] = {0.0};
G4double G4BGGPionElasticXS::theGlauberFacPiMinus[93] = {0.0};
G4double G4BGGPionElasticXS::theCoulombFacPiMinus[93] = {0.0};
G4int    G4BGGPionElasticXS::theA[93] = {0};

#ifdef G4MULTITHREADED
G4Mutex G4BGGPionElasticXS::pionElasticXSMutex = G4MUTEX_INITIALIZER;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BGGPionElasticXS::G4BGGPionElasticXS(const G4ParticleDefinition* p) 
 : G4VCrossSectionDataSet("BarashenkovGlauberGribov") 
{
  verboseLevel = 0;
  fGlauberEnergy = 91.*GeV;
  fLowEnergy = 20.*MeV;
  fLowestEnergy = 1.*MeV;
  SetMinKinEnergy(0.0);
  SetMaxKinEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );

  fPion = nullptr;
  fGlauber = nullptr;
  fHadron  = nullptr;

  fG4pow   = G4Pow::GetInstance();

  theProton= G4Proton::Proton();
  thePiPlus= G4PionPlus::PionPlus();
  isPiplus = (p == thePiPlus);
  isMaster = false;
  SetForAllAtomsAndEnergies(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4BGGPionElasticXS::~G4BGGPionElasticXS()
{
  delete fHadron;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool 
G4BGGPionElasticXS::IsElementApplicable(const G4DynamicParticle*, G4int,
                                        const G4Material*)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4BGGPionElasticXS::IsIsoApplicable(const G4DynamicParticle*, 
                                           G4int Z, G4int, 
                                           const G4Element*, const G4Material*)
{
  return (1 == Z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4BGGPionElasticXS::GetElementCrossSection(const G4DynamicParticle* dp,
                                           G4int ZZ, const G4Material*)
{
  // this method should be called only for Z > 1
  G4double cross = 0.0;
  G4double ekin = std::max(dp->GetKineticEnergy(), fLowestEnergy);
  G4int Z = std::min(ZZ, 92);
  if(1 == Z) {
    cross = 1.0115*GetIsoCrossSection(dp,1,1);
  } else {
    if(ekin <= fLowEnergy) {
      cross = (isPiplus) ? theCoulombFacPiPlus[Z]*CoulombFactorPiPlus(ekin, Z)
        : theCoulombFacPiMinus[Z]*FactorPiMinus(ekin);
    } else if(ekin > fGlauberEnergy) {
      cross = (isPiplus) ? theGlauberFacPiPlus[Z] : theGlauberFacPiMinus[Z];
      cross *= fGlauber->GetElasticGlauberGribov(dp, Z, theA[Z]);
    } else {
      cross = fPion->GetElasticCrossSection(dp, Z, theA[Z]);
    }
  }
  if(verboseLevel > 1) {
    G4cout << "G4BGGPionElasticXS::GetElementCrossSection  for "
           << dp->GetDefinition()->GetParticleName()
           << "  Ekin(GeV)= " << dp->GetKineticEnergy()
           << " in nucleus Z= " << Z << "  A= " << theA[Z]
           << " XS(b)= " << cross/barn 
           << G4endl;
  }
  return cross;
}

G4double
G4BGGPionElasticXS::GetIsoCrossSection(const G4DynamicParticle* dp, 
                                       G4int Z, G4int A, 
                                       const G4Isotope*,
                                       const G4Element*,
                                       const G4Material*)
{
  // this method should be called only for Z = 1
  fHadron->HadronNucleonXscNS(dp->GetDefinition(), theProton, 
                              dp->GetKineticEnergy());
  G4double cross = A*fHadron->GetElasticHadronNucleonXsc();

  if(verboseLevel > 1) {
    G4cout << "G4BGGPionElasticXS::GetIsoCrossSection  for "
           << dp->GetDefinition()->GetParticleName()
           << "  Ekin(GeV)= " << dp->GetKineticEnergy()
           << " in nucleus Z= " << Z << "  A= " << A
           << " XS(b)= " << cross/barn 
           << G4endl;
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4BGGPionElasticXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(fPion) { return; }
  if(verboseLevel > 1) {
    G4cout << "G4BGGPionElasticXS::BuildPhysicsTable for " 
           << p.GetParticleName() << G4endl;
  } 
  if(&p == G4PionPlus::PionPlus() || &p == G4PionMinus::PionMinus()) {
    isPiplus = (&p == G4PionPlus::PionPlus());
  } else {
    G4ExceptionDescription ed;
    ed << "This BGG cross section is applicable only to pions and not to " 
       << p.GetParticleName() << G4endl; 
    G4Exception("G4BGGPionElasticXS::BuildPhysicsTable", "had001", 
                FatalException, ed);
    return;
  }

  fPion    = new G4UPiNuclearCrossSection();
  fGlauber = new G4ComponentGGHadronNucleusXsc();
  fHadron  = new G4HadronNucleonXsc();

  fPion->BuildPhysicsTable(p);

  if(0 == theA[0]) { 
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&pionElasticXSMutex);
    if(0 == theA[0]) { 
#endif
      isMaster = true;
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&pionElasticXSMutex);
#endif
  } else {
    return;
  }

  if(isMaster && 0 == theA[0]) {

    theA[0] = theA[1] = 1;
    G4ThreeVector mom(0.0,0.0,1.0);
    G4DynamicParticle dp(thePiPlus, mom, fGlauberEnergy);

    G4NistManager* nist = G4NistManager::Instance();

    G4double csup, csdn;
    for(G4int iz=2; iz<93; ++iz) {

      G4int A = G4lrint(nist->GetAtomicMassAmu(iz));
      theA[iz] = A;

      csup = fGlauber->GetElasticGlauberGribov(&dp, iz, A);
      csdn = fPion->GetElasticCrossSection(&dp, iz, A);
      theGlauberFacPiPlus[iz] = csdn/csup;
    }

    dp.SetDefinition(G4PionMinus::PionMinus());
    for(G4int iz=2; iz<93; ++iz) {
      csup = fGlauber->GetElasticGlauberGribov(&dp, iz, theA[iz]);
      csdn = fPion->GetElasticCrossSection(&dp, iz, theA[iz]);
      theGlauberFacPiMinus[iz] = csdn/csup;

      if(verboseLevel > 0) {
        G4cout << "Z= " << iz <<  "  A= " << theA[iz] 
               << " factorPiPlus= " << theGlauberFacPiPlus[iz] 
               << " factorPiMinus= " << theGlauberFacPiMinus[iz] 
               << G4endl; 
      }
    }
    theCoulombFacPiPlus[1] = 1.0;
    theCoulombFacPiMinus[1]= 1.0;
    dp.SetKineticEnergy(fLowEnergy);
    dp.SetDefinition(thePiPlus);
    for(G4int iz=2; iz<93; ++iz) {
      theCoulombFacPiPlus[iz] = fPion->GetElasticCrossSection(&dp, iz, theA[iz])
        /CoulombFactorPiPlus(fLowEnergy, iz);
    }
    dp.SetDefinition(G4PionMinus::PionMinus());
    for(G4int iz=2; iz<93; ++iz) {
      theCoulombFacPiMinus[iz] = fPion->GetElasticCrossSection(&dp, iz, theA[iz])
        /FactorPiMinus(fLowEnergy);

      if(verboseLevel > 0) {
        G4cout << "Z= " << iz <<  "  A= " << theA[iz] 
               << " CoulombFactorPiPlus= " << theCoulombFacPiPlus[iz] 
               << " CoulombFactorPiMinus= " << theCoulombFacPiMinus[iz] 
               << G4endl; 
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BGGPionElasticXS::CoulombFactorPiPlus(G4double kinEnergy, G4int Z)
{
  return (kinEnergy > 0.0) ? 
    G4NuclearRadii::CoulombFactor(Z, theA[Z], thePiPlus, kinEnergy) : 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4BGGPionElasticXS::FactorPiMinus(G4double kinEnergy)
{
  return 1.0/std::sqrt(kinEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
G4BGGPionElasticXS::CrossSectionDescription(std::ostream& outFile) const 
{
  outFile << "The Barashenkov-Glauber-Gribov cross section handles elastic\n"
          << "scattering of pions from nuclei at all energies. The\n"
          << "Barashenkov parameterization is used below 91 GeV and the\n"
          << "Glauber-Gribov parameterization is used above 91 GeV.\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
