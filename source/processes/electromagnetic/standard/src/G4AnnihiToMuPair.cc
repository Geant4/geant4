//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4AnnihiToMuPair.cc,v 1.2 2003-02-04 11:08:42 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//         ------------ G4AnnihiToMuPair physics process ------
//         by H.Burkhardt, S. Kelner and R. Kokoulin, November 2002
// -----------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......//
//
// 04.02.03 : cosmetic simplifications (mma)
// 27.01.03 : first implementation (hbu)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4AnnihiToMuPair.hh"

#include "G4ios.hh"
#include "Randomize.hh"

#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Material.hh"
#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// constructor

G4AnnihiToMuPair::G4AnnihiToMuPair(const G4String& processName)
  : G4VDiscreteProcess (processName)
{
 //e+ Energy threshold
 const G4double Mu_massc2 = G4MuonPlus::MuonPlus()->GetPDGMass();
 LowestEnergyLimit  = 2*Mu_massc2*Mu_massc2/electron_mass_c2 - electron_mass_c2;
 
 //modele ok up to 1000 TeV due to neglected Z-interference
 HighestEnergyLimit = 1000*TeV;
 
 CrossSecFactor = 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4AnnihiToMuPair::~G4AnnihiToMuPair() // (empty) destructor
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4AnnihiToMuPair::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Positron::Positron() );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AnnihiToMuPair::BuildPhysicsTable(const G4ParticleDefinition&)
// Build cross section and mean free path tables
//here no tables, just calling PrintInfoDefinition
{
   PrintInfoDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AnnihiToMuPair::SetCrossSecFactor(G4double fac)
// Set the factor to artificially increase the cross section
{ 
  CrossSecFactor = fac;
  G4cout << "The cross section for AnnihiToMuPair is artificially "
         << "increased by the CrossSecFactor=" << CrossSecFactor << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4AnnihiToMuPair::ComputeCrossSectionPerAtom(G4double Epos, G4double Z)
// Calculates the microscopic cross section in GEANT4 internal units.
// It gives a good description from threshold to 1000 GeV
{
  static const G4double Mmuon = G4MuonPlus::MuonPlus()->GetPDGMass();
  static const G4double Rmuon = elm_coupling/Mmuon; //classical particle radius
  static const G4double Sig0  = pi*Rmuon*Rmuon/3.;  //constant in crossSection

  G4double CrossSection = 0.;
  if (Epos < LowestEnergyLimit) return CrossSection;
   
  G4double xi = LowestEnergyLimit/Epos;
  G4double SigmaEl = Sig0*xi*(1.+xi/2.)*sqrt(1.-xi); // per electron
  CrossSection = SigmaEl*Z;         // number of electrons per atom
  CrossSection *= CrossSecFactor;   //increase the CrossSection by  (default 1)
  return CrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4AnnihiToMuPair::GetMeanFreePath(const G4Track& aTrack,
                                           G4double, G4ForceCondition*)

// returns the positron mean free path in GEANT4 internal units

{
  const G4DynamicParticle* aDynamicPositron = aTrack.GetDynamicParticle();
  G4double PositronEnergy = aDynamicPositron->GetKineticEnergy()
                                              +electron_mass_c2;
  G4Material* aMaterial = aTrack.GetMaterial();
  const G4ElementVector* theElementVector = aMaterial->GetElementVector();
  const G4double* NbOfAtomsPerVolume = aMaterial->GetVecNbOfAtomsPerVolume();

  G4double SIGMA = 0 ;

  for ( size_t i=0 ; i < aMaterial->GetNumberOfElements() ; i++ )
  {
    G4double AtomicZ = (*theElementVector)[i]->GetZ();
    SIGMA += NbOfAtomsPerVolume[i] *
      ComputeCrossSectionPerAtom(PositronEnergy,AtomicZ);
  }
  return SIGMA > DBL_MIN ? 1./SIGMA : DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* G4AnnihiToMuPair::PostStepDoIt(const G4Track& aTrack,
                                                  const G4Step&  aStep)
//
// generation of e+e- -> mu+mu-
//
{

  aParticleChange.Initialize(aTrack);
  static const G4double Mele=electron_mass_c2;
  static const G4double Mmuon=G4MuonPlus::MuonPlus()->GetPDGMass();

  // current Positron energy and direction, return if energy too low
  const G4DynamicParticle *aDynamicPositron = aTrack.GetDynamicParticle();
  G4double Epos = aDynamicPositron->GetKineticEnergy()+Mele;

 if (Epos < LowestEnergyLimit)
  { G4cout 
        << "error in G4AnnihiToMuPair::PostStepDoIt called with energy below"
	   " threshold Epos= "
	<< Epos << G4endl;       // shoud never happen
	G4Exception(10);
  }

  if (Epos < LowestEnergyLimit)
     return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);

  G4ParticleMomentum PositronDirection = 
                                       aDynamicPositron->GetMomentumDirection();
  G4double xi = LowestEnergyLimit/Epos; // xi is always less than 1,
                                        // goes to 0 at high Epos

  // generate cost
  //
  G4double cost;
  do cost = 2.*G4UniformRand()-1.;
  while (2.*G4UniformRand() > 1.+xi+cost*cost*(1.-xi) ); 
                                                       //1+cost**2 at high Epos
  G4double sint = sqrt(1.-cost*cost);

  // generate phi
  //
  G4double phi=2.*pi*G4UniformRand();

  G4double Ecm   = sqrt(0.5*Mele*(Epos+Mele));
  G4double Pcm   = sqrt(Ecm*Ecm-Mmuon*Mmuon);
  G4double beta  = sqrt((Epos-Mele)/(Epos+Mele));
  G4double gamma = Ecm/Mele;                    // =sqrt((Epos+Mele)/(2.*Mele));
  G4double Pt    = Pcm*sint;
  
  // energy and momentum of the muons in the Lab
  //
  G4double EmuPlus   = gamma*(     Ecm+cost*beta*Pcm);
  G4double EmuMinus  = gamma*(     Ecm-cost*beta*Pcm);
  G4double PmuPlusZ  = gamma*(beta*Ecm+cost*     Pcm);
  G4double PmuMinusZ = gamma*(beta*Ecm-cost*     Pcm);
  G4double PmuPlusX  = Pt*cos(phi);
  G4double PmuPlusY  = Pt*sin(phi);
  G4double PmuMinusX =-Pt*cos(phi);
  G4double PmuMinusY =-Pt*sin(phi);
  // absolute momenta
  G4double PmuPlus  = sqrt(Pt*Pt+PmuPlusZ *PmuPlusZ );
  G4double PmuMinus = sqrt(Pt*Pt+PmuMinusZ*PmuMinusZ);

  // mu+ mu- directions for Positron in z-direction
  //
  G4ThreeVector
    MuPlusDirection ( PmuPlusX/PmuPlus, PmuPlusY/PmuPlus,  PmuPlusZ/PmuPlus  );
  G4ThreeVector
    MuMinusDirection(PmuMinusX/PmuMinus,PmuMinusY/PmuMinus,PmuMinusZ/PmuMinus);

  // rotate to actual Positron direction
  //
  MuPlusDirection.rotateUz(PositronDirection);
  MuMinusDirection.rotateUz(PositronDirection);

  aParticleChange.SetNumberOfSecondaries(2);
  // create G4DynamicParticle object for the particle1
  G4DynamicParticle* aParticle1= new G4DynamicParticle(
                         G4MuonPlus::MuonPlus(),MuPlusDirection,EmuPlus-Mmuon);
  aParticleChange.AddSecondary(aParticle1);
  // create G4DynamicParticle object for the particle2
  G4DynamicParticle* aParticle2= new G4DynamicParticle(
                     G4MuonMinus::MuonMinus(),MuMinusDirection,EmuMinus-Mmuon);
  aParticleChange.AddSecondary(aParticle2);

  // Kill the incident positron 
  //
  aParticleChange.SetMomentumChange( 0., 0., 0. );
  aParticleChange.SetEnergyChange(0.); 
  aParticleChange.SetStatusChange(fStopAndKill);

  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AnnihiToMuPair::PrintInfoDefinition()
{
  G4String comments ="e+e->mu+mu- annihilation, atomic e- at rest.\n";
  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "        threshold at " << LowestEnergyLimit/GeV << " GeV"
         << " good description up to "
         << HighestEnergyLimit/TeV << " TeV for all Z." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
