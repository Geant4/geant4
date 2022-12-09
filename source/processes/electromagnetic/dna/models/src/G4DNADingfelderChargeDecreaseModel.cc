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

#include "G4DNADingfelderChargeDecreaseModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4DNAChemistryManager.hh"
#include "G4Log.hh"
#include "G4Pow.hh"
#include "G4Alpha.hh"


static G4Pow * gpow = G4Pow::GetInstance();

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNADingfelderChargeDecreaseModel::G4DNADingfelderChargeDecreaseModel(const G4ParticleDefinition*,
                                                                       const G4String& nam) :
G4VEmModel(nam), isInitialised(false)
{
  numberOfPartialCrossSections[0] = 0;
  numberOfPartialCrossSections[1] = 0;
  numberOfPartialCrossSections[2] = 0;

  verboseLevel = 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if (verboseLevel > 0)
  {
    G4cout << "Dingfelder charge decrease model is constructed " << G4endl;
  }
  // Selection of stationary mode

  statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNADingfelderChargeDecreaseModel::Initialise(const G4ParticleDefinition* particle,
                                                    const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
  {
    G4cout << "Calling G4DNADingfelderChargeDecreaseModel::Initialise()"
        << G4endl;
  }

  // Energy limits

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  protonDef = G4Proton::ProtonDefinition();
  alphaPlusPlusDef = G4Alpha::Alpha();
  alphaPlusDef = instance->GetIon("alpha+");
  hydrogenDef = instance->GetIon("hydrogen");
  heliumDef = instance->GetIon("helium");

  G4String proton;
  G4String alphaPlusPlus;
  G4String alphaPlus;

  // Limits

  proton = protonDef->GetParticleName();
  lowEnergyLimit[proton] = 100. * eV;
  highEnergyLimit[proton] = 100. * MeV;

  alphaPlusPlus = alphaPlusPlusDef->GetParticleName();
  lowEnergyLimit[alphaPlusPlus] = 1. * keV;
  highEnergyLimit[alphaPlusPlus] = 400. * MeV;

  alphaPlus = alphaPlusDef->GetParticleName();
  lowEnergyLimit[alphaPlus] = 1. * keV;
  highEnergyLimit[alphaPlus] = 400. * MeV;

  //

  if (particle==protonDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[proton]);
    SetHighEnergyLimit(highEnergyLimit[proton]);
  }

  if (particle==alphaPlusPlusDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[alphaPlusPlus]);
    SetHighEnergyLimit(highEnergyLimit[alphaPlusPlus]);
  }

  if (particle==alphaPlusDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[alphaPlus]);
    SetHighEnergyLimit(highEnergyLimit[alphaPlus]);
  }

  // Final state

  //PROTON
  f0[0][0]=1.;
  a0[0][0]=-0.180;
  a1[0][0]=-3.600;
  b0[0][0]=-18.22;
  b1[0][0]=-1.997;
  c0[0][0]=0.215;
  d0[0][0]=3.550;
  x0[0][0]=3.450;
  x1[0][0]=5.251;

  numberOfPartialCrossSections[0] = 1;

  //ALPHA++
  f0[0][1]=1.; a0[0][1]=0.95;
  a1[0][1]=-2.75;
  b0[0][1]=-23.00;
  c0[0][1]=0.215;
  d0[0][1]=2.95;
  x0[0][1]=3.50;

  f0[1][1]=1.;
  a0[1][1]=0.95;
  a1[1][1]=-2.75;
  b0[1][1]=-23.73;
  c0[1][1]=0.250;
  d0[1][1]=3.55;
  x0[1][1]=3.72;

  x1[0][1]=-1.;
  b1[0][1]=-1.;

  x1[1][1]=-1.;
  b1[1][1]=-1.;

  numberOfPartialCrossSections[1] = 2;

  // ALPHA+
  f0[0][2]=1.;
  a0[0][2]=0.65;
  a1[0][2]=-2.75;
  b0[0][2]=-21.81;
  c0[0][2]=0.232;
  d0[0][2]=2.95;
  x0[0][2]=3.53;

  x1[0][2]=-1.;
  b1[0][2]=-1.;

  numberOfPartialCrossSections[2] = 1;

  //

  if( verboseLevel>0 )
  {
    G4cout << "Dingfelder charge decrease model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / keV << " keV - "
    << HighEnergyLimit() / MeV << " MeV for "
    << particle->GetParticleName()
    << G4endl;
  }

  // Initialize water density pointer
  fpMolWaterDensity = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

  if (isInitialised)
  { return;}
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNADingfelderChargeDecreaseModel::CrossSectionPerVolume(const G4Material* material,
                                                                   const G4ParticleDefinition* particleDefinition,
                                                                   G4double k,
                                                                   G4double,
                                                                   G4double)
{
  if (verboseLevel > 3)
  {
    G4cout
        << "Calling CrossSectionPerVolume() of G4DNADingfelderChargeDecreaseModel"
        << G4endl;
  }

  // Calculate total cross section for model
  if (
      particleDefinition != protonDef
      &&
      particleDefinition != alphaPlusPlusDef
      &&
      particleDefinition != alphaPlusDef   
  )

  return 0;

  G4double lowLim = 0;
  G4double highLim = 0;
  G4double crossSection = 0.;

  G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];

  const G4String& particleName = particleDefinition->GetParticleName();

  std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
  pos1 = lowEnergyLimit.find(particleName);

  if (pos1 != lowEnergyLimit.end())
  {
    lowLim = pos1->second;
  }

  std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
  pos2 = highEnergyLimit.find(particleName);

  if (pos2 != highEnergyLimit.end())
  {
    highLim = pos2->second;
  }

  if (k >= lowLim && k <= highLim)
  {
    crossSection = Sum(k,particleDefinition);
  }

  if (verboseLevel > 2)
  {
    G4cout << "_______________________________________" << G4endl;
    G4cout << "G4DNADingfelderChargeDecreaeModel" << G4endl;
    G4cout << "Kinetic energy(eV)=" << k/eV << "particle :" << particleName << G4endl;
    G4cout << "Cross section per water molecule (cm^2)=" << crossSection/cm/cm << G4endl;
    G4cout << "Cross section per water molecule (cm^-1)=" << crossSection*
    waterDensity/(1./cm) << G4endl;
    // material->GetAtomicNumDensityVector()[1]/(1./cm) << G4endl;
  }

  return crossSection*waterDensity;
  //return crossSection*material->GetAtomicNumDensityVector()[1];

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNADingfelderChargeDecreaseModel::SampleSecondaries(std::vector<
                                                               G4DynamicParticle*>* fvect,
                                                           const G4MaterialCutsCouple* /*couple*/,
                                                           const G4DynamicParticle* aDynamicParticle,
                                                           G4double,
                                                           G4double)
{
  if (verboseLevel > 3)
  {
    G4cout
        << "Calling SampleSecondaries() of G4DNADingfelderChargeDecreaseModel"
        << G4endl;
  }

  G4double inK = aDynamicParticle->GetKineticEnergy();

  G4ParticleDefinition* definition = aDynamicParticle->GetDefinition();

  G4double particleMass = definition->GetPDGMass();

  G4int finalStateIndex = RandomSelect(inK,definition);

  G4int n = NumberOfFinalStates(definition, finalStateIndex);
  G4double waterBindingEnergy = WaterBindingEnergyConstant(definition, finalStateIndex);
  G4double outgoingParticleBindingEnergy = OutgoingParticleBindingEnergyConstant(definition, finalStateIndex);

  G4double outK = 0.;
  
  if (definition==G4Proton::Proton())
  {
   if (!statCode) outK = inK - n*(inK*electron_mass_c2/proton_mass_c2) 
                             - waterBindingEnergy + outgoingParticleBindingEnergy;
   else outK = inK;
  }
  
  else
  {
   if (!statCode) outK = inK - n*(inK*electron_mass_c2/particleMass) 
                             - waterBindingEnergy + outgoingParticleBindingEnergy;
   else outK = inK;
  }
  
  if (outK<0)
  {
    G4Exception("G4DNADingfelderChargeDecreaseModel::SampleSecondaries","em0004",
        FatalException,"Final kinetic energy is negative.");
  }

  fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);

  if (!statCode) fParticleChangeForGamma->ProposeLocalEnergyDeposit(waterBindingEnergy);
  
  else
  
  {
   if (definition==G4Proton::Proton()) 
     fParticleChangeForGamma->ProposeLocalEnergyDeposit(n*(inK*electron_mass_c2/proton_mass_c2) 
     + waterBindingEnergy - outgoingParticleBindingEnergy);
   else
     fParticleChangeForGamma->ProposeLocalEnergyDeposit(n*(inK*electron_mass_c2/particleMass) 
     + waterBindingEnergy - outgoingParticleBindingEnergy);
  }

  G4DynamicParticle* dp = new G4DynamicParticle (OutgoingParticleDefinition(definition, finalStateIndex),
      aDynamicParticle->GetMomentumDirection(),
      outK);
  fvect->push_back(dp);

  const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
      G4DNAChemistryManager::Instance()->CreateWaterMolecule(eIonizedMolecule,
          1,
          theIncomingTrack);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNADingfelderChargeDecreaseModel::NumberOfFinalStates(G4ParticleDefinition* particleDefinition,
                                                              G4int finalStateIndex)

{
  if (particleDefinition == G4Proton::Proton())
    return 1;

  if (particleDefinition == alphaPlusPlusDef)
  {
    if (finalStateIndex == 0)
      return 1;
    return 2;
  }

  if (particleDefinition == alphaPlusDef)
    return 1;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleDefinition* G4DNADingfelderChargeDecreaseModel::OutgoingParticleDefinition(G4ParticleDefinition* particleDefinition,
                                                                                     G4int finalStateIndex)
{
  if (particleDefinition == G4Proton::Proton())
    return hydrogenDef;

  if (particleDefinition == alphaPlusPlusDef)
  {
    if (finalStateIndex == 0)
      return alphaPlusDef;
    return heliumDef;
  }

  if (particleDefinition == alphaPlusDef)
    return heliumDef;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNADingfelderChargeDecreaseModel::WaterBindingEnergyConstant(G4ParticleDefinition* particleDefinition,
                                                                        G4int finalStateIndex)
{
  // Ionization energy of first water shell
  // Rad. Phys. Chem. 59 p.267 by Dingf. et al.
  // W + 10.79 eV -> W+ + e-

  if (particleDefinition == G4Proton::Proton())
    return 10.79 * eV;

  if (particleDefinition == alphaPlusPlusDef)
  {
    // Binding energy for    W+ -> W++ + e-    10.79 eV
    // Binding energy for    W  -> W+  + e-    10.79 eV
    //
    // Ionization energy of first water shell
    // Rad. Phys. Chem. 59 p.267 by Dingf. et al.

    if (finalStateIndex == 0)
      return 10.79 * eV;

    return 10.79 * 2 * eV;
  }

  if (particleDefinition == alphaPlusDef)
  {
    // Binding energy for    W+ -> W++ + e-    10.79 eV
    // Binding energy for    W  -> W+  + e-    10.79 eV
    //
    // Ionization energy of first water shell
    // Rad. Phys. Chem. 59 p.267 by Dingf. et al.

    return 10.79 * eV;
  }

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNADingfelderChargeDecreaseModel::OutgoingParticleBindingEnergyConstant(G4ParticleDefinition* particleDefinition,
                                                                                   G4int finalStateIndex)
{
  if (particleDefinition == G4Proton::Proton())
    return 13.6 * eV;

  if (particleDefinition == alphaPlusPlusDef)
  {
    // Binding energy for    He+ -> He++ + e-    54.509 eV
    // Binding energy for    He  -> He+  + e-    24.587 eV

    if (finalStateIndex == 0)
      return 54.509 * eV;

    return (54.509 + 24.587) * eV;
  }

  if (particleDefinition == alphaPlusDef)
  {
    // Binding energy for    He+ -> He++ + e-    54.509 eV
    // Binding energy for    He  -> He+  + e-    24.587 eV

    return 24.587 * eV;
  }

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNADingfelderChargeDecreaseModel::PartialCrossSection(G4double k,
                                                                 G4int index,
                                                                 const G4ParticleDefinition* particleDefinition)
{
  G4int particleTypeIndex = 0;

  if (particleDefinition == protonDef)
    particleTypeIndex = 0;

  if (particleDefinition == alphaPlusPlusDef)
    particleTypeIndex = 1;

  if (particleDefinition == alphaPlusDef)
    particleTypeIndex = 2;

  //
  // sigma(T) = f0 10 ^ y(log10(T/eV))
  //
  //         /  a0 x + b0                    x < x0
  //         |
  // y(x) = <   a0 x + b0 - c0 (x - x0)^d0   x0 <= x < x1
  //         |
  //         \  a1 x + b1                    x >= x1
  //
  //
  // f0, a0, a1, b0, b1, c0, d0, x0, x1 are parameters that change for protons and helium (0, +, ++)
  //
  // f0 has been added to the code in order to manage partial (shell-dependent) cross sections (if no shell dependence is present. f0=1. Sum of f0 over the considered shells should give 1)
  //
  // From Rad. Phys. and Chem. 59 (2000) 255-275, M. Dingfelder et al.
  // Inelastic-collision cross sections of liquid water for interactions of energetic proton
  //

  if (x1[index][particleTypeIndex] < x0[index][particleTypeIndex])
  {
    //
    // if x1 < x0 means that x1 and b1 will be calculated with the following formula (this piece of code is run on all alphas and not on protons)
    //
    // x1 = x0 + ((a0 - a1)/(c0 * d0)) ^ (1 / (d0 - 1))
    //
    // b1 = (a0 - a1) * x1 + b0 - c0 * (x1 - x0) ^ d0
    //

    x1[index][particleTypeIndex] = x0[index][particleTypeIndex]
        + gpow->powA((a0[index][particleTypeIndex] - a1[index][particleTypeIndex])
                       / (c0[index][particleTypeIndex]
                           * d0[index][particleTypeIndex]),
                   1. / (d0[index][particleTypeIndex] - 1.));
    b1[index][particleTypeIndex] = (a0[index][particleTypeIndex]
        - a1[index][particleTypeIndex]) * x1[index][particleTypeIndex]
        + b0[index][particleTypeIndex]
        - c0[index][particleTypeIndex]
            * gpow->powA(x1[index][particleTypeIndex]
                           - x0[index][particleTypeIndex],
                       d0[index][particleTypeIndex]);
  }

  G4double x(G4Log(k / eV)/gpow->logZ(10));
  G4double y;

  if (x < x0[index][particleTypeIndex])
    y = a0[index][particleTypeIndex] * x + b0[index][particleTypeIndex];
  else if (x < x1[index][particleTypeIndex])
    y = a0[index][particleTypeIndex] * x + b0[index][particleTypeIndex]
        - c0[index][particleTypeIndex]
            * gpow->powA(x - x0[index][particleTypeIndex],
                       d0[index][particleTypeIndex]);
  else
    y = a1[index][particleTypeIndex] * x + b1[index][particleTypeIndex];

  return f0[index][particleTypeIndex] * gpow->powA(10., y) * m * m;

}

G4int G4DNADingfelderChargeDecreaseModel::RandomSelect(G4double k,
                                                       const G4ParticleDefinition* particleDefinition)
{
  G4int particleTypeIndex = 0;

  if (particleDefinition == protonDef)
    particleTypeIndex = 0;

  if (particleDefinition == alphaPlusPlusDef)
    particleTypeIndex = 1;

  if (particleDefinition == alphaPlusDef)
    particleTypeIndex = 2;

  const G4int n = numberOfPartialCrossSections[particleTypeIndex];
  G4double* values(new G4double[n]);
  G4double value(0);
  G4int i = n;

  while (i > 0)
  {
    i--;
    values[i] = PartialCrossSection(k, i, particleDefinition);
    value += values[i];
  }

  value *= G4UniformRand();

  i = n;
  while (i > 0)
  {
    i--;

    if (values[i] > value)
      break;

    value -= values[i];
  }

  delete[] values;

  return i;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNADingfelderChargeDecreaseModel::Sum(G4double k,
                                                 const G4ParticleDefinition* particleDefinition)
{
  G4int particleTypeIndex = 0;

  if (particleDefinition == protonDef)
    particleTypeIndex = 0;

  if (particleDefinition == alphaPlusPlusDef)
    particleTypeIndex = 1;

  if (particleDefinition == alphaPlusDef)
    particleTypeIndex = 2;

  G4double totalCrossSection = 0.;

  for (G4int i = 0; i < numberOfPartialCrossSections[particleTypeIndex]; i++)
  {
    totalCrossSection += PartialCrossSection(k, i, particleDefinition);
  }

  return totalCrossSection;
}

