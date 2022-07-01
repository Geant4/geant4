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

#include "G4DNADingfelderChargeIncreaseModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4Log.hh"
#include "G4Pow.hh"
#include "G4Alpha.hh"

static G4Pow * gpow = G4Pow::GetInstance();

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNADingfelderChargeIncreaseModel::G4DNADingfelderChargeIncreaseModel(const G4ParticleDefinition*,
                                                                       const G4String& nam) :
G4VEmModel(nam), isInitialised(false)
{
  fpMolWaterDensity = 0;

  numberOfPartialCrossSections[0] = 0;
  numberOfPartialCrossSections[1] = 0;

  verboseLevel = 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if (verboseLevel > 0)
  {
    G4cout << "Dingfelder charge increase model is constructed " << G4endl;
  }
  fParticleChangeForGamma = 0;

  // Selection of stationary mode

  statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNADingfelderChargeIncreaseModel::~G4DNADingfelderChargeIncreaseModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNADingfelderChargeIncreaseModel::Initialise(const G4ParticleDefinition* particle,
                                                    const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
  {
    G4cout << "Calling G4DNADingfelderChargeIncreaseModel::Initialise()"
        << G4endl;
  }

  // Energy limits

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  hydrogenDef = instance->GetIon("hydrogen");
  alphaPlusPlusDef = G4Alpha::Alpha();
  alphaPlusDef = instance->GetIon("alpha+");
  heliumDef = instance->GetIon("helium");

  G4String hydrogen;
  G4String alphaPlus;
  G4String helium;

  // Limits

  hydrogen = hydrogenDef->GetParticleName();
  lowEnergyLimit[hydrogen] = 100. * eV;
  highEnergyLimit[hydrogen] = 100. * MeV;

  alphaPlus = alphaPlusDef->GetParticleName();
  lowEnergyLimit[alphaPlus] = 1. * keV;
  highEnergyLimit[alphaPlus] = 400. * MeV;

  helium = heliumDef->GetParticleName();
  lowEnergyLimit[helium] = 1. * keV;
  highEnergyLimit[helium] = 400. * MeV;

  //

  if (particle==hydrogenDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[hydrogen]);
    SetHighEnergyLimit(highEnergyLimit[hydrogen]);
  }

  if (particle==alphaPlusDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[alphaPlus]);
    SetHighEnergyLimit(highEnergyLimit[alphaPlus]);
  }

  if (particle==heliumDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[helium]);
    SetHighEnergyLimit(highEnergyLimit[helium]);
  }

  // Final state

  //ALPHA+

  f0[0][0]=1.;
  a0[0][0]=2.25;
  a1[0][0]=-0.75;
  b0[0][0]=-32.10;
  c0[0][0]=0.600;
  d0[0][0]=2.40;
  x0[0][0]=4.60;

  x1[0][0]=-1.;
  b1[0][0]=-1.;

  numberOfPartialCrossSections[0]=1;

  //HELIUM

  f0[0][1]=1.;
  a0[0][1]=2.25;
  a1[0][1]=-0.75;
  b0[0][1]=-30.93;
  c0[0][1]=0.590;
  d0[0][1]=2.35;
  x0[0][1]=4.29;

  f0[1][1]=1.;
  a0[1][1]=2.25;
  a1[1][1]=-0.75;
  b0[1][1]=-32.61;
  c0[1][1]=0.435;
  d0[1][1]=2.70;
  x0[1][1]=4.45;

  x1[0][1]=-1.;
  b1[0][1]=-1.;

  x1[1][1]=-1.;
  b1[1][1]=-1.;

  numberOfPartialCrossSections[1]=2;

  //

  if( verboseLevel>0 )
  {
    G4cout << "Dingfelder charge increase model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / keV << " keV - "
    << HighEnergyLimit() / MeV << " MeV for "
    << particle->GetParticleName()
    << G4endl;
  }

  // Initialize water density pointer
  fpMolWaterDensity = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

  if (isInitialised) return;

  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNADingfelderChargeIncreaseModel::CrossSectionPerVolume(const G4Material* material,
                                                                   const G4ParticleDefinition* particleDefinition,
                                                                   G4double k,
                                                                   G4double,
                                                                   G4double)
{
  if (verboseLevel > 3)
  {
    G4cout
        << "Calling CrossSectionPerVolume() of G4DNADingfelderChargeIncreaseModel"
        << G4endl;
  }

  // Calculate total cross section for model

  if (
      particleDefinition != hydrogenDef
      &&
      particleDefinition != alphaPlusDef
      &&
      particleDefinition != heliumDef
  )

  return 0;

  G4double lowLim = 0;
  G4double highLim = 0;
  G4double totalCrossSection = 0.;

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
    //HYDROGEN
    if (particleDefinition == hydrogenDef)
    {
      const G4double aa = 2.835;
      const G4double bb = 0.310;
      const G4double cc = 2.100;
      const G4double dd = 0.760;
      const G4double fac = 1.0e-18;
      const G4double rr = 13.606 * eV;

      G4double t = k / (proton_mass_c2/electron_mass_c2);
      G4double x = t / rr;
      G4double temp = 4.0 * pi * Bohr_radius/nm * Bohr_radius/nm * fac;
      G4double sigmal = temp * cc * (gpow->powA(x,dd));
      G4double sigmah = temp * (aa * G4Log(1.0 + x) + bb) / x;
      totalCrossSection = 1.0/(1.0/sigmal + 1.0/sigmah) *m*m;
    }
    else
    {
      totalCrossSection = Sum(k,particleDefinition);
    }
  }

  if (verboseLevel > 2)
  {
    G4cout << "__________________________________" << G4endl;
    G4cout << "G4DNADingfelderChargeIncreaseModel - XS INFO START" << G4endl;
    G4cout << "Kinetic energy(eV)=" << k/eV << " particle : " << particleName << G4endl;
    G4cout << "Cross section per water molecule (cm^2)=" << totalCrossSection/cm/cm << G4endl;
    G4cout << "Cross section per water molecule (cm^-1)=" << totalCrossSection*waterDensity/(1./cm) << G4endl;
    //  G4cout << " - Cross section per water molecule (cm^-1)=" 
    //  << sigma*material->GetAtomicNumDensityVector()[1]/(1./cm) << G4endl;
    G4cout << "G4DNADingfelderChargeIncreaseModel - XS INFO END" << G4endl;
  }

  return totalCrossSection*waterDensity;
  // return totalCrossSection*material->GetAtomicNumDensityVector()[1];

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNADingfelderChargeIncreaseModel::SampleSecondaries(std::vector<
                                                               G4DynamicParticle*>* fvect,
                                                           const G4MaterialCutsCouple* /*couple*/,
                                                           const G4DynamicParticle* aDynamicParticle,
                                                           G4double,
                                                           G4double)
{
  if (verboseLevel > 3)
  {
    G4cout
        << "Calling SampleSecondaries() of G4DNADingfelderChargeIncreaseModel"
        << G4endl;
  }
  
  if (!statCode) fParticleChangeForGamma->ProposeLocalEnergyDeposit(0.);
  
  G4ParticleDefinition* definition = aDynamicParticle->GetDefinition();

  G4double particleMass = definition->GetPDGMass();

  G4double inK = aDynamicParticle->GetKineticEnergy();

  G4int finalStateIndex = RandomSelect(inK,definition);

  G4int n = NumberOfFinalStates(definition,finalStateIndex);

  G4double outK = 0.;
  
  if (!statCode) outK = inK - IncomingParticleBindingEnergyConstant(definition,finalStateIndex);

  else outK = inK;
  
  if (statCode) 
    fParticleChangeForGamma->
      ProposeLocalEnergyDeposit(IncomingParticleBindingEnergyConstant(definition,finalStateIndex));

  fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);

  G4double electronK;
  if (definition == hydrogenDef) electronK = inK*electron_mass_c2/proton_mass_c2;
  else electronK = inK*electron_mass_c2/(particleMass);

  if (outK<0)
  {
    G4Exception("G4DNADingfelderChargeIncreaseModel::SampleSecondaries","em0004",
        FatalException,"Final kinetic energy is negative.");
  }

  G4DynamicParticle* dp = new G4DynamicParticle(OutgoingParticleDefinition(definition,finalStateIndex),
      aDynamicParticle->GetMomentumDirection(),
      outK);

  fvect->push_back(dp);

  n = n - 1;

  while (n>0)
  {
    n--;
    fvect->push_back(new G4DynamicParticle
        (G4Electron::Electron(), aDynamicParticle->GetMomentumDirection(), electronK) );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNADingfelderChargeIncreaseModel::NumberOfFinalStates(G4ParticleDefinition* particleDefinition,
                                                              G4int finalStateIndex)

{

  if (particleDefinition == hydrogenDef)
    return 2;

  if (particleDefinition == alphaPlusDef)
    return 2;

  if (particleDefinition == heliumDef)
  {
    if (finalStateIndex == 0)
      return 2;
    return 3;
  }

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleDefinition* G4DNADingfelderChargeIncreaseModel::OutgoingParticleDefinition(G4ParticleDefinition* particleDefinition,
                                                                                     G4int finalStateIndex)
{

  if (particleDefinition == hydrogenDef)
    return G4Proton::Proton();

  if (particleDefinition == alphaPlusDef)
    return alphaPlusPlusDef;

  if (particleDefinition == heliumDef)
  {
    if (finalStateIndex == 0)
      return alphaPlusDef;
    return alphaPlusPlusDef;
  }

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNADingfelderChargeIncreaseModel::IncomingParticleBindingEnergyConstant(G4ParticleDefinition* particleDefinition,
                                                                                   G4int finalStateIndex)
{

  if (particleDefinition == hydrogenDef)
    return 13.6 * eV;

  if (particleDefinition == alphaPlusDef)
  {
    // Binding energy for    He+ -> He++ + e-    54.509 eV
    // Binding energy for    He  -> He+  + e-    24.587 eV
    return 54.509 * eV;
  }

  if (particleDefinition == heliumDef)
  {
    // Binding energy for    He+ -> He++ + e-    54.509 eV
    // Binding energy for    He  -> He+  + e-    24.587 eV

    if (finalStateIndex == 0)
      return 24.587 * eV;
    return (54.509 + 24.587) * eV;
  }

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNADingfelderChargeIncreaseModel::PartialCrossSection(G4double k,
                                                                 G4int index,
                                                                 const G4ParticleDefinition* particleDefinition)
{
  G4int particleTypeIndex = 0;

  if (particleDefinition == alphaPlusDef)
    particleTypeIndex = 0;

  if (particleDefinition == heliumDef)
    particleTypeIndex = 1;

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
  // f0 has been added to the code in order to manage partial (shell-dependent) cross sections 
  // (if no shell dependence is present. f0=1. Sum of f0 over the considered shells should give 1)
  //
  // From Rad. Phys. and Chem. 59 (2000) 255-275, M. Dingfelder et al.
  // Inelastic-collision cross sections of liquid water for interactions of energetic proton
  //

  if (x1[index][particleTypeIndex] < x0[index][particleTypeIndex])
  {
    //
    // if x1 < x0 means that x1 and b1 will be calculated with the following formula 
    // (this piece of code is run on all alphas and not on protons)
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNADingfelderChargeIncreaseModel::RandomSelect(G4double k,
                                                       const G4ParticleDefinition* particleDefinition)
{
  G4int particleTypeIndex = 0;

  if (particleDefinition == hydrogenDef)
    return 0;

  if (particleDefinition == alphaPlusDef)
    particleTypeIndex = 0;

  if (particleDefinition == heliumDef)
    particleTypeIndex = 1;

  const G4int n = numberOfPartialCrossSections[particleTypeIndex];
  G4double* values(new G4double[n]);
  G4double value = 0;
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

G4double G4DNADingfelderChargeIncreaseModel::Sum(G4double k,
                                                 const G4ParticleDefinition* particleDefinition)
{
  G4int particleTypeIndex = 0;

  if (particleDefinition == alphaPlusDef)
    particleTypeIndex = 0;

  if (particleDefinition == heliumDef)
    particleTypeIndex = 1;

  G4double totalCrossSection = 0.;

  for (G4int i = 0; i < numberOfPartialCrossSections[particleTypeIndex]; i++)
  {
    totalCrossSection += PartialCrossSection(k, i, particleDefinition);
  }

  return totalCrossSection;
}
