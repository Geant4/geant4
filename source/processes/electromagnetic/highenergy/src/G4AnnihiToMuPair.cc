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
//
//         ------------ G4AnnihiToMuPair physics process ------
//         by H.Burkhardt, S. Kelner and R. Kokoulin, November 2002
// -----------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......//
//
// 27.01.03 : first implementation (hbu)
// 04.02.03 : cosmetic simplifications (mma)
// 25.10.04 : migrade to new interfaces of ParticleChange (vi)
// 28.02.18 : cross section now including SSS threshold factor
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4AnnihiToMuPair.hh"

#include "G4ios.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4TauPlus.hh"
#include "G4TauMinus.hh"
#include "G4Material.hh"
#include "G4Step.hh"
#include "G4LossTableManager.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4AnnihiToMuPair::G4AnnihiToMuPair(const G4String& processName,
    G4ProcessType type):G4VDiscreteProcess (processName, type)
{
  //e+ Energy threshold
  if(processName == "AnnihiToTauPair") {
    SetProcessSubType(fAnnihilationToTauTau);
    part1 = G4TauPlus::TauPlus();
    part2 = G4TauMinus::TauMinus();
    fInfo = "e+e->tau+tau-";
  } else {
    SetProcessSubType(fAnnihilationToMuMu);
    part1 = G4MuonPlus::MuonPlus();
    part2 = G4MuonMinus::MuonMinus();
  }
  fMass = part1->GetPDGMass();
  fLowEnergyLimit = 
    2.*fMass*fMass/CLHEP::electron_mass_c2 - CLHEP::electron_mass_c2;
 
  //model is ok up to 1000 TeV due to neglected Z-interference
  fHighEnergyLimit = 1000.*TeV;
 
  fCurrentSigma = 0.0;
  fCrossSecFactor = 1.;
  fManager = G4LossTableManager::Instance();
  fManager->Register(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4AnnihiToMuPair::~G4AnnihiToMuPair() // (empty) destructor
{ 
  fManager->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4AnnihiToMuPair::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Positron::Positron() );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AnnihiToMuPair::BuildPhysicsTable(const G4ParticleDefinition&)
{
  PrintInfoDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AnnihiToMuPair::SetCrossSecFactor(G4double fac)
// Set the factor to artificially increase the cross section
{ 
  fCrossSecFactor = fac;
  //G4cout << "The cross section for AnnihiToMuPair is artificially "
  //       << "increased by the CrossSecFactor=" << fCrossSecFactor << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4AnnihiToMuPair::ComputeCrossSectionPerElectron(const G4double e)
// Calculates the microscopic cross section in GEANT4 internal units.
// It gives a good description from threshold to 1000 GeV
{
  G4double rmuon = CLHEP::elm_coupling/fMass; //classical particle radius
  G4double sig0 = CLHEP::pi*rmuon*rmuon/3.;   //constant in crossSection
  const G4double pial = CLHEP::pi*CLHEP::fine_structure_const; // pi * alphaQED

  if (e <= fLowEnergyLimit) return 0.0;
   
  const G4double xi = fLowEnergyLimit/e;
  const G4double piaxi = pial * std::sqrt(xi);
  G4double sigma = sig0 * xi * (1. + xi*0.5);
  //G4cout << "### xi= " << xi << " piaxi=" << piaxi << G4endl;

  // argument of the exponent below 0.1 or above 10
  // Sigma per electron * number of electrons per atom
  if(xi <= 1.0 - 100*piaxi*piaxi) {
    sigma *= std::sqrt(1.0 - xi);
  } else if( xi >= 1.0 - 0.01*piaxi*piaxi) {
    sigma *= piaxi;
  } else {
    sigma *= piaxi/(1. - G4Exp( -piaxi/std::sqrt(1-xi) ));
  }
  //G4cout << "### sigma= " << sigma << G4endl;
  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4AnnihiToMuPair::ComputeCrossSectionPerAtom(const G4double energy,
                                                      const G4double Z)
{
  return ComputeCrossSectionPerElectron(energy)*Z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4AnnihiToMuPair::CrossSectionPerVolume(G4double energy, 
						 const G4Material* aMaterial)
{
  return ComputeCrossSectionPerElectron(energy)*aMaterial->GetTotNbOfElectPerVolume();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4AnnihiToMuPair::GetMeanFreePath(const G4Track& aTrack,
                                           G4double, G4ForceCondition*)
// returns the positron mean free path in GEANT4 internal units
{
  const G4DynamicParticle* aDynamicPositron = aTrack.GetDynamicParticle();
  G4double energy = aDynamicPositron->GetTotalEnergy();
  const G4Material* aMaterial = aTrack.GetMaterial();

  // cross section before step
  fCurrentSigma = CrossSectionPerVolume(energy, aMaterial);

  // increase the CrossSection by CrossSecFactor (default 1)
  return (fCurrentSigma > 0.0) ? 1.0/(fCurrentSigma*fCrossSecFactor) : DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* G4AnnihiToMuPair::PostStepDoIt(const G4Track& aTrack,
                                                  const G4Step&  aStep)
//
// generation of e+e- -> mu+mu-
//
{
  aParticleChange.Initialize(aTrack);

  // current Positron energy and direction, return if energy too low
  const G4DynamicParticle *aDynamicPositron = aTrack.GetDynamicParticle();
  const G4double Mele = CLHEP::electron_mass_c2;
  G4double Epos = aDynamicPositron->GetTotalEnergy();
  G4double xs = CrossSectionPerVolume(Epos, aTrack.GetMaterial());

  // test of cross section
  if(xs > 0.0 && fCurrentSigma*G4UniformRand() > xs) 
    {
      return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
    }

  const G4ThreeVector PosiDirection = aDynamicPositron->GetMomentumDirection();
  G4double xi = fLowEnergyLimit/Epos; // xi is always less than 1,
                                      // goes to 0 at high Epos

  // generate cost; probability function 1+cost**2 at high Epos
  //
  G4double cost;
  do { cost = 2.*G4UniformRand()-1.; }
  // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  while (2.*G4UniformRand() > 1.+xi+cost*cost*(1.-xi) ); 
  G4double sint = sqrt(1.-cost*cost);

  // generate phi
  //
  G4double phi = 2.*CLHEP::pi*G4UniformRand();

  G4double Ecm   = std::sqrt(0.5*Mele*(Epos+Mele));
  G4double Pcm   = std::sqrt(Ecm*Ecm - fMass*fMass);
  G4double beta  = std::sqrt((Epos-Mele)/(Epos+Mele));
  G4double gamma = Ecm/Mele;
  G4double Pt    = Pcm*sint;
  
  // energy and momentum of the muons in the Lab
  //
  G4double EmuPlus   = gamma*(Ecm + cost*beta*Pcm);
  G4double EmuMinus  = gamma*(Ecm - cost*beta*Pcm);
  G4double PmuPlusZ  = gamma*(beta*Ecm + cost*Pcm);
  G4double PmuMinusZ = gamma*(beta*Ecm - cost*Pcm);
  G4double PmuPlusX  = Pt*std::cos(phi);
  G4double PmuPlusY  = Pt*std::sin(phi);
  G4double PmuMinusX =-PmuPlusX;
  G4double PmuMinusY =-PmuPlusY;
  // absolute momenta
  G4double PmuPlus  = std::sqrt(Pt*Pt+PmuPlusZ *PmuPlusZ );
  G4double PmuMinus = std::sqrt(Pt*Pt+PmuMinusZ*PmuMinusZ);

  // mu+ mu- directions for Positron in z-direction
  //
  G4ThreeVector
    MuPlusDirection(PmuPlusX/PmuPlus, PmuPlusY/PmuPlus, PmuPlusZ/PmuPlus);
  G4ThreeVector
    MuMinusDirection(PmuMinusX/PmuMinus,PmuMinusY/PmuMinus,PmuMinusZ/PmuMinus);

  // rotate to actual Positron direction
  //
  MuPlusDirection.rotateUz(PosiDirection);
  MuMinusDirection.rotateUz(PosiDirection);

  aParticleChange.SetNumberOfSecondaries(2);

  // create G4DynamicParticle object for the particle1
  G4DynamicParticle* aParticle1 = 
    new G4DynamicParticle(part1, MuPlusDirection, EmuPlus-fMass);
  aParticleChange.AddSecondary(aParticle1);
  // create G4DynamicParticle object for the particle2
  G4DynamicParticle* aParticle2 =
    new G4DynamicParticle(part2, MuMinusDirection, EmuMinus-fMass);
  aParticleChange.AddSecondary(aParticle2);

  // Kill the incident positron 
  //
  aParticleChange.ProposeEnergy(0.); 
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AnnihiToMuPair::PrintInfoDefinition()
{
  G4String comments = fInfo + " annihilation, atomic e- at rest, SubType=";
  G4cout << G4endl << GetProcessName() << ":  " << comments 
	 << GetProcessSubType() << G4endl;
  G4cout << "        threshold at " << fLowEnergyLimit/CLHEP::GeV << " GeV"
         << " good description up to "
         << fHighEnergyLimit/CLHEP::TeV << " TeV for all Z." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
