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
// $Id: G4DNABornExcitationModel1.cc 90054 2015-05-11 18:47:32Z matkara $
//

#include "G4DNABornExcitationModel1.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNABornExcitationModel1::G4DNABornExcitationModel1(const G4ParticleDefinition*,
                                                     const G4String& nam) :
    G4VEmModel(nam), isInitialised(false), fTableData(0)
{
  fpMolWaterDensity = 0;
  fHighEnergy = 0;
  fLowEnergy = 0;
  fParticleDefinition = 0;

  verboseLevel = 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if (verboseLevel > 0)
  {
    G4cout << "Born excitation model is constructed " << G4endl;
  }
  fParticleChangeForGamma = 0;

  // Selection of stationary mode

  statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNABornExcitationModel1::~G4DNABornExcitationModel1()
{
  // Cross section
  if (fTableData)
    delete fTableData;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNABornExcitationModel1::Initialise(const G4ParticleDefinition* particle,
                                           const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
  {
    G4cout << "Calling G4DNABornExcitationModel1::Initialise()" << G4endl;
  }

  if(fParticleDefinition != 0 && fParticleDefinition != particle)
  {
    G4Exception("G4DNABornExcitationModel1::Initialise","em0001",
        FatalException,"Model already initialized for another particle type.");
  }

  fParticleDefinition = particle;

  if(particle->GetParticleName() == "e-")
  {
    fTableFile = "dna/sigma_excitation_e_born";
    fLowEnergy = 9*eV;
    fHighEnergy = 1*MeV;
  }
  else if(particle->GetParticleName() == "proton")
  {
    fTableFile = "dna/sigma_excitation_p_born";
    fLowEnergy = 500. * keV;
    fHighEnergy = 100. * MeV;
  }

  SetLowEnergyLimit(fLowEnergy);
  SetHighEnergyLimit(fHighEnergy);

  G4double scaleFactor = (1.e-22 / 3.343) * m*m;
  fTableData = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
  fTableData->LoadData(fTableFile);

  if( verboseLevel>0 )
  {
    G4cout << "Born excitation model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / eV << " eV - "
    << HighEnergyLimit() / keV << " keV for "
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

G4double G4DNABornExcitationModel1::CrossSectionPerVolume(const G4Material* material,
                                                          const G4ParticleDefinition* particleDefinition,
                                                          G4double ekin,
                                                          G4double,
                                                          G4double)
{
  if (verboseLevel > 3)
  {
    G4cout << "Calling CrossSectionPerVolume() of G4DNABornExcitationModel1"
        << G4endl;
  }

  if(particleDefinition != fParticleDefinition) return 0;

  // Calculate total cross section for model

  G4double sigma=0;

  G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];

  if(waterDensity!= 0.0)
  {
    if (ekin >= fLowEnergy && ekin < fHighEnergy)
    {
      sigma = fTableData->FindValue(ekin);
    }

    if (verboseLevel > 2)
    {
      G4cout << "__________________________________" << G4endl;
      G4cout << "G4DNABornExcitationModel1 - XS INFO START" << G4endl;
      G4cout << "Kinetic energy(eV)=" << ekin/eV << " particle : " << particleDefinition->GetParticleName() << G4endl;
      G4cout << "Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
      G4cout << "Cross section per water molecule (cm^-1)=" << sigma*waterDensity/(1./cm) << G4endl;
      G4cout << "G4DNABornExcitationModel1 - XS INFO END" << G4endl;
    }
  } // if (waterMaterial)

  return sigma*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNABornExcitationModel1::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
                                                  const G4MaterialCutsCouple* /*couple*/,
                                                  const G4DynamicParticle* aDynamicParticle,
                                                  G4double,
                                                  G4double)
{

  if (verboseLevel > 3)
  {
    G4cout << "Calling SampleSecondaries() of G4DNABornExcitationModel1"
        << G4endl;
  }

  G4double k = aDynamicParticle->GetKineticEnergy();

  G4int level = RandomSelect(k);
  G4double excitationEnergy = waterStructure.ExcitationEnergy(level);
  G4double newEnergy = k - excitationEnergy;

  if (newEnergy > 0)
  {
    fParticleChangeForGamma->ProposeMomentumDirection(aDynamicParticle->GetMomentumDirection());
    
    if (!statCode) fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
    else fParticleChangeForGamma->SetProposedKineticEnergy(k);
    
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
  }

  const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
  G4DNAChemistryManager::Instance()->CreateWaterMolecule(eExcitedMolecule,
      level,
      theIncomingTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNABornExcitationModel1::GetPartialCrossSection(const G4Material*,
                                                           G4int level,
                                                           const G4ParticleDefinition* particle,
                                                           G4double kineticEnergy)
{
  if (fParticleDefinition != particle)
  {
    G4Exception("G4DNABornExcitationModel1::GetPartialCrossSection",
                "bornParticleType",
                FatalException,
                "Model initialized for another particle type.");
  }

  return fTableData->GetComponent(level)->FindValue(kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNABornExcitationModel1::RandomSelect(G4double k)
{
  G4int level = 0;

  G4double* valuesBuffer = new G4double[fTableData->NumberOfComponents()];
  const size_t n(fTableData->NumberOfComponents());
  size_t i(n);
  G4double value = 0.;

  while (i > 0)
  {
    i--;
    valuesBuffer[i] = fTableData->GetComponent(i)->FindValue(k);
    value += valuesBuffer[i];
  }

  value *= G4UniformRand();
  i = n;

  while (i > 0)
  {
    i--;

    if (valuesBuffer[i] > value)
    {
      delete[] valuesBuffer;
      return i;
    }
    value -= valuesBuffer[i];
  }

  if (valuesBuffer)
    delete[] valuesBuffer;

  return level;
}

