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
// Reference:
//    A.D. Dominguez-Munoz, M.I. Gallardo, M.C. Bordage,
//    Z. Francis, S. Incerti, M.A. Cortes-Giraldo,
//    Radiat. Phys. Chem. 199 (2022) 110363.
//
// Class authors:
//    A.D. Dominguez-Munoz
//    M.A. Cortes-Giraldo (miancortes -at- us.es)
//
// Class creation: 2022-03-03
//
//

#include "G4DNARPWBAExcitationModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include <map>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARPWBAExcitationModel::G4DNARPWBAExcitationModel(
  const G4ParticleDefinition*, const G4String& nam)
  : G4VEmModel(nam)
{
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods
  if(verboseLevel > 0)
  {
    G4cout << "RPWBA excitation model is constructed " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARPWBAExcitationModel::~G4DNARPWBAExcitationModel() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARPWBAExcitationModel::Initialise(const G4ParticleDefinition* particle,
                                           const G4DataVector& /*cuts*/)
{
  if(isInitialised)
  {
    return;
  }
  if(verboseLevel > 3)
  {
    G4cout << "Calling G4DNARPWBAExcitationModel::Initialise()" << G4endl;
  }

  if(fParticleDefinition != nullptr && fParticleDefinition != particle)
  {
    G4Exception("G4DNARPWBAExcitationModel::Initialise", "em0001",
                FatalException,
                "Model already initialized for another particle type.");
  }

  fTableFile  = "dna/sigma_excitation_p_RPWBA";
  fLowEnergy  = 100. * MeV;
  fHighEnergy = 300. * MeV;

  if(LowEnergyLimit() < fLowEnergy || HighEnergyLimit() > fHighEnergy)
  {
    G4ExceptionDescription ed;
    ed << "Model is applicable from "<<fLowEnergy<<" to "<<fHighEnergy;
    G4Exception("G4DNARPWBAExcitationModel::Initialise", "em0004",
      FatalException, ed);
  }

  G4double scaleFactor = 1 * cm * cm;
  fTableData = make_unique<G4DNACrossSectionDataSet>(new G4LogLogInterpolation,
                                                     eV, scaleFactor);
  fTableData->LoadData(fTableFile);

  if(verboseLevel > 0)
  {
    G4cout << "RPWBA excitation model is initialized " << G4endl
           << "Energy range: " << LowEnergyLimit() / eV << " eV - "
           << HighEnergyLimit() / keV << " keV for "
           << particle->GetParticleName() << G4endl;
  }

  // Initialize water density pointer
  if(G4Material::GetMaterial("G4_WATER") != nullptr){
    fpMolWaterDensity =
      G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(
        G4Material::GetMaterial("G4_WATER"));
  }else{
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "G4_WATER does not exist :";
    G4Exception("G4DNARPWBAIonisationModel::Initialise", "em00020",
                FatalException, exceptionDescription);
  }
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised           = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNARPWBAExcitationModel::CrossSectionPerVolume(
  const G4Material* material, const G4ParticleDefinition* particleDefinition,
  G4double ekin, G4double, G4double)
{
  if(verboseLevel > 3)
  {
    G4cout << "Calling CrossSectionPerVolume() of G4DNARPWBAExcitationModel"
           << G4endl;
  }

  if(fTableData == nullptr)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "No cross section data ";
    G4Exception("G4DNARPWBAIonisationModel::CrossSectionPerVolume", "em00120",
                FatalException, exceptionDescription);
  }

  if(particleDefinition != fParticleDefinition)
    return 0;

  // Calculate total cross section for model

  G4double sigma = 0;

  G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];

  if(ekin >= fLowEnergy && ekin <= fHighEnergy)
  {
    sigma = fTableData->FindValue(ekin);
  }

  if(verboseLevel > 2)
  {
    G4cout << "__________________________________" << G4endl;
    G4cout << "G4DNARPWBAExcitationModel - XS INFO START" << G4endl;
    G4cout << "Kinetic energy(eV)=" << ekin / eV
           << " particle : " << particleDefinition->GetParticleName() << G4endl;
    G4cout << "Cross section per water molecule (cm^2)=" << sigma / cm / cm
           << G4endl;
    G4cout << "Cross section per water molecule (cm^-1)="
           << sigma * waterDensity / (1. / cm) << G4endl;
    G4cout << "G4DNARPWBAExcitationModel - XS INFO END" << G4endl;
  }

  return sigma * waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARPWBAExcitationModel::SampleSecondaries(
  std::vector<G4DynamicParticle*>* /*fvect*/,
  const G4MaterialCutsCouple* /*couple*/,
  const G4DynamicParticle* aDynamicParticle, G4double, G4double)
{
  if(verboseLevel > 3)
  {
    G4cout << "Calling SampleSecondaries() of G4DNARPWBAExcitationModel"
           << G4endl;
  }

  G4double k = aDynamicParticle->GetKineticEnergy();

  G4int level               = RandomSelect(k);
  G4double excitationEnergy = waterStructure.ExcitationEnergy(level);
  G4double newEnergy        = k - excitationEnergy;

  if(newEnergy > 0)
  {
    fParticleChangeForGamma->ProposeMomentumDirection(
      aDynamicParticle->GetMomentumDirection());

    if(!statCode){
      fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
    }
    else{
      fParticleChangeForGamma->SetProposedKineticEnergy(k);
    }
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
  }

  const G4Track* theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
  G4DNAChemistryManager::Instance()->CreateWaterMolecule(
    eExcitedMolecule, level, theIncomingTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARPWBAExcitationModel::GetPartialCrossSection(
  const G4Material*, G4int level, const G4ParticleDefinition* particle,
  G4double kineticEnergy)
{
  if(fParticleDefinition != particle)
  {
    G4Exception("G4DNARPWBAExcitationModel::GetPartialCrossSection",
                "RPWBAParticleType", FatalException,
                "Model initialized for another particle type.");
  }

  return fTableData->GetComponent(level)->FindValue(kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNARPWBAExcitationModel::RandomSelect(G4double k)
{
  G4int level = 0;

  G4double* valuesBuffer = new G4double[fTableData->NumberOfComponents()];
  const G4int n = (G4int)fTableData->NumberOfComponents();
  G4int i(n);
  G4double value = 0.;

  while(i > 0)
  {
    --i;
    valuesBuffer[i] = fTableData->GetComponent(i)->FindValue(k);
    value += valuesBuffer[i];
  }

  value *= G4UniformRand();
  i = n;

  while(i > 0)
  {
    --i;

    if(valuesBuffer[i] > value)
    {
      delete[] valuesBuffer;
      return i;
    }
    value -= valuesBuffer[i];
  }
  delete[] valuesBuffer;

  return level;
}
