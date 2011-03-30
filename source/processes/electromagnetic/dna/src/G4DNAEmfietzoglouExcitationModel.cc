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
// $Id: G4DNAEmfietzoglouExcitationModel.cc,v 1.10 2010-06-08 21:50:00 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4DNAEmfietzoglouExcitationModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAEmfietzoglouExcitationModel::G4DNAEmfietzoglouExcitationModel(const G4ParticleDefinition*,
                                             const G4String& nam)
:G4VEmModel(nam),isInitialised(false)
{

  lowEnergyLimit = 8.23 * eV; 
  highEnergyLimit = 10 * MeV;
  SetLowEnergyLimit(lowEnergyLimit);
  SetHighEnergyLimit(highEnergyLimit);

  nLevels = waterExcitation.NumberOfLevels();

  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods
  
  if (verboseLevel > 3)

  if( verboseLevel>0 ) 
  { 
    G4cout << "Emfietzoglou Excitation model is constructed " << G4endl
           << "Energy range: "
           << lowEnergyLimit / eV << " eV - "
           << highEnergyLimit / MeV << " MeV"
           << G4endl;
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAEmfietzoglouExcitationModel::~G4DNAEmfietzoglouExcitationModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAEmfietzoglouExcitationModel::Initialise(const G4ParticleDefinition* /*particle*/,
                                       const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
    G4cout << "Calling G4DNAEmfietzoglouExcitationModel::Initialise()" << G4endl;

  // Energy limits
  
  if (LowEnergyLimit() < lowEnergyLimit)
  {
    G4cout << "G4DNAEmfietzoglouExcitationModel: low energy limit increased from " << 
	LowEnergyLimit()/eV << " eV to " << lowEnergyLimit/eV << " eV" << G4endl;
    SetLowEnergyLimit(lowEnergyLimit);
    }

  if (HighEnergyLimit() > highEnergyLimit)
  {
    G4cout << "G4DNAEmfietzoglouExcitationModel: high energy limit decreased from " << 
        HighEnergyLimit()/MeV << " MeV to " << highEnergyLimit/MeV << " MeV" << G4endl;
    SetHighEnergyLimit(highEnergyLimit);
  }

  //
  if( verboseLevel>0 ) 
  { 
    G4cout << "Emfietzoglou Excitation model is initialized " << G4endl
           << "Energy range: "
           << LowEnergyLimit() / eV << " eV - "
           << HighEnergyLimit() / MeV << " MeV"
           << G4endl;
  }

  if(!isInitialised) 
  {
    isInitialised = true;
  
    if(pParticleChange)
      fParticleChangeForGamma = reinterpret_cast<G4ParticleChangeForGamma*>(pParticleChange);
    else
      fParticleChangeForGamma = new G4ParticleChangeForGamma();
  }    

  // InitialiseElementSelectors(particle,cuts);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAEmfietzoglouExcitationModel::CrossSectionPerVolume(const G4Material* material,
					   const G4ParticleDefinition* particleDefinition,
					   G4double ekin,
					   G4double,
					   G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of G4DNAEmfietzoglouExcitationModel" << G4endl;

 // Calculate total cross section for model

 G4double sigma=0;
 
 if (material->GetName() == "G4_WATER")
 {

  if (particleDefinition == G4Electron::ElectronDefinition())
  {
    if (ekin >= lowEnergyLimit && ekin < highEnergyLimit)
    {
      sigma = Sum(ekin);
    }
  }
  
  if (verboseLevel > 3)
  {
    G4cout << "---> Kinetic energy(eV)=" << ekin/eV << G4endl;
    G4cout << " - Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
    G4cout << " - Cross section per water molecule (cm^-1)=" << sigma*material->GetAtomicNumDensityVector()[1]/(1./cm) << G4endl;
  } 

 } 
       
 return sigma*material->GetAtomicNumDensityVector()[1];		   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAEmfietzoglouExcitationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
					      const G4MaterialCutsCouple* /*couple*/,
					      const G4DynamicParticle* aDynamicElectron,
					      G4double,
					      G4double)
{

  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4DNAEmfietzoglouExcitationModel" << G4endl;

  G4double electronEnergy0 = aDynamicElectron->GetKineticEnergy();
  
  G4int level = RandomSelect(electronEnergy0);

  G4double excitationEnergy = waterExcitation.ExcitationEnergy(level);
  G4double newEnergy = electronEnergy0 - excitationEnergy;
  
  if (electronEnergy0 < highEnergyLimit)
  {
    if (newEnergy >= lowEnergyLimit)
    {
      fParticleChangeForGamma->ProposeMomentumDirection(aDynamicElectron->GetMomentumDirection());
      fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
    }
 
    else   
    {
      fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(electronEnergy0);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNAEmfietzoglouExcitationModel::PartialCrossSection(G4double t, G4int level)
{
  //                 Aj                        T
  // Sigma(T) = ------------- (Bj /  T) ln(Cj ---) [1 - Bj / T]^Pj
  //             2 pi alpha0                   R
  //
  // Sigma is the macroscopic cross section = N sigma, where N = number of target particles per unit volume
  // and sigma is the microscopic cross section
  // T      is the incoming electron kinetic energy
  // alpha0 is the Bohr Radius (Bohr_radius)
  // Aj, Bj, Cj & Pj are parameters that can be found in Emfietzoglou's papers
  //
  // From Phys. Med. Biol. 48 (2003) 2355-2371, D.Emfietzoglou,
  // Monte Carlo Simulation of the energy loss of low energy electrons in liquid Water
  //
  // Scaling for macroscopic cross section: number of water moleculs per unit volume
  // const G4double sigma0 = (10. / 3.343e22) * cm2;

  const G4double density = 3.34192e+19 * mm3;

  const G4double aj[]={0.0205, 0.0209, 0.0130, 0.0026, 0.0025};
  const G4double cj[]={4.9801, 3.3850, 2.8095, 1.9242, 3.4624};
  const G4double pj[]={0.4757, 0.3483, 0.4443, 0.3429, 0.4379};
  const G4double r = 13.6 * eV;
  
  G4double sigma = 0.;
  
  G4double exc = waterExcitation.ExcitationEnergy(level);
  
  if (t >= exc)
  {
      G4double excitationSigma = ( aj[level] / (2.*pi*Bohr_radius)) 
	* (exc / t) 
	* std::log(cj[level]*(t/r)) 
	* std::pow((1.- (exc/t)), pj[level]);
      sigma = excitationSigma / density;
  }

  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNAEmfietzoglouExcitationModel::RandomSelect(G4double k)
{
  G4int i = nLevels;
  G4double value = 0.;
  std::deque<double> values;
  
  while (i > 0)
  {
    i--;
    G4double partial = PartialCrossSection(k,i);
    values.push_front(partial);
    value += partial;
  }

  value *= G4UniformRand();
    
  i = nLevels;

  while (i > 0)
  {
    i--;
    if (values[i] > value) return i;
    value -= values[i];
  }
    
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNAEmfietzoglouExcitationModel::Sum(G4double k)
{
  G4double totalCrossSection = 0.;

  for (G4int i=0; i<nLevels; i++)
  {
    totalCrossSection += PartialCrossSection(k,i);
  }
  return totalCrossSection;
}

