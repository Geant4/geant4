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
// Geant4 Class file
//
// File name:     G4PolarizedIonisation
//
// Author:        A.Schaelicke on base of Vladimir Ivanchenko code

#include "G4PolarizedIonisation.hh"

#include "G4Electron.hh"
#include "G4EmParameters.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4PolarizationHelper.hh"
#include "G4PolarizationManager.hh"
#include "G4PolarizedIonisationModel.hh"
#include "G4Positron.hh"
#include "G4ProductionCutsTable.hh"
#include "G4StokesVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4UniversalFluctuation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4PolarizedIonisation::G4PolarizedIonisation(const G4String& name)
  : G4VEnergyLossProcess(name)
  , fAsymmetryTable(nullptr)
  , fTransverseAsymmetryTable(nullptr)
  , fIsElectron(true)
  , fIsInitialised(false)
{
  verboseLevel = 0;
  SetProcessSubType(fIonisation);
  SetSecondaryParticle(G4Electron::Electron());
  fFlucModel = nullptr;
  fEmModel   = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4PolarizedIonisation::~G4PolarizedIonisation() { CleanTables(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedIonisation::ProcessDescription(std::ostream& out) const
{
  out << "Polarized version of G4eIonisation.\n";

  G4VEnergyLossProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedIonisation::CleanTables()
{
  if(fAsymmetryTable)
  {
    fAsymmetryTable->clearAndDestroy();
    delete fAsymmetryTable;
    fAsymmetryTable = nullptr;
  }
  if(fTransverseAsymmetryTable)
  {
    fTransverseAsymmetryTable->clearAndDestroy();
    delete fTransverseAsymmetryTable;
    fTransverseAsymmetryTable = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4PolarizedIonisation::MinPrimaryEnergy(const G4ParticleDefinition*,
                                                 const G4Material*,
                                                 G4double cut)
{
  G4double x = cut;
  if(fIsElectron)
  {
    x += cut;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4bool G4PolarizedIonisation::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Electron::Electron() || &p == G4Positron::Positron());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedIonisation::InitialiseEnergyLossProcess(
  const G4ParticleDefinition* part, const G4ParticleDefinition*)
{
  if(!fIsInitialised)
  {
    if(part == G4Positron::Positron())
    {
      fIsElectron = false;
    }

    if(!FluctModel())
    {
      SetFluctModel(new G4UniversalFluctuation());
    }
    fFlucModel = FluctModel();

    fEmModel = new G4PolarizedIonisationModel();
    SetEmModel(fEmModel);
    G4EmParameters* param = G4EmParameters::Instance();
    fEmModel->SetLowEnergyLimit(param->MinKinEnergy());
    fEmModel->SetHighEnergyLimit(param->MaxKinEnergy());
    AddEmModel(1, fEmModel, fFlucModel);

    fIsInitialised = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4PolarizedIonisation::GetMeanFreePath(const G4Track& track,
                                                G4double step,
                                                G4ForceCondition* cond)
{
  // *** get unploarised mean free path from lambda table ***
  G4double mfp = G4VEnergyLossProcess::GetMeanFreePath(track, step, cond);
  if(fAsymmetryTable && fTransverseAsymmetryTable && mfp < DBL_MAX)
  {
    mfp *= ComputeSaturationFactor(track);
  }
  if(verboseLevel >= 2)
  {
    G4cout << "G4PolarizedIonisation::MeanFreePath:  " << mfp / mm << " mm "
           << G4endl;
  }
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4PolarizedIonisation::PostStepGetPhysicalInteractionLength(
  const G4Track& track, G4double step, G4ForceCondition* cond)
{
  // save previous values
  G4double nLength = theNumberOfInteractionLengthLeft;
  G4double iLength = currentInteractionLength;

  // *** get unpolarised mean free path from lambda table ***
  // this changes theNumberOfInteractionLengthLeft and currentInteractionLength
  G4double x = G4VEnergyLossProcess::PostStepGetPhysicalInteractionLength(
    track, step, cond);
  G4double x0      = x;
  G4double satFact = 1.;

  // *** add corrections on polarisation ***
  if(fAsymmetryTable && fTransverseAsymmetryTable && x < DBL_MAX)
  {
    satFact            = ComputeSaturationFactor(track);
    G4double curLength = currentInteractionLength * satFact;
    G4double prvLength = iLength * satFact;
    if(nLength > 0.0)
    {
      theNumberOfInteractionLengthLeft =
        std::max(nLength - step / prvLength, 0.0);
    }
    x = theNumberOfInteractionLengthLeft * curLength;
  }
  if(verboseLevel >= 2)
  {
    G4cout << "G4PolarizedIonisation::PostStepGPIL: " << std::setprecision(8)
           << x / mm << " mm;" << G4endl
           << "                   unpolarized value: " << std::setprecision(8)
           << x0 / mm << " mm." << G4endl;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4PolarizedIonisation::ComputeSaturationFactor(const G4Track& track)
{
  const G4Material* aMaterial = track.GetMaterial();
  G4VPhysicalVolume* aPVolume = track.GetVolume();
  G4LogicalVolume* aLVolume   = aPVolume->GetLogicalVolume();

  G4PolarizationManager* polarizationManager =
    G4PolarizationManager::GetInstance();

  const G4bool volumeIsPolarized = polarizationManager->IsPolarized(aLVolume);
  G4StokesVector volPolarization =
    polarizationManager->GetVolumePolarization(aLVolume);

  G4double factor = 1.0;

  if(volumeIsPolarized && !volPolarization.IsZero())
  {
    // *** get asymmetry, if target is polarized ***
    const G4DynamicParticle* aDynamicPart = track.GetDynamicParticle();
    const G4double energy                 = aDynamicPart->GetKineticEnergy();
    const G4StokesVector polarization = G4StokesVector(track.GetPolarization());
    const G4ParticleMomentum direction0 = aDynamicPart->GetMomentumDirection();

    if(verboseLevel >= 2)
    {
      G4cout << "G4PolarizedIonisation::ComputeSaturationFactor: " << G4endl;
      G4cout << " Energy(MeV)  " << energy / MeV << G4endl;
      G4cout << " Direction    " << direction0 << G4endl;
      G4cout << " Polarization " << polarization << G4endl;
      G4cout << " MaterialPol. " << volPolarization << G4endl;
      G4cout << " Phys. Volume " << aPVolume->GetName() << G4endl;
      G4cout << " Log. Volume  " << aLVolume->GetName() << G4endl;
      G4cout << " Material     " << aMaterial << G4endl;
    }

    std::size_t midx               = CurrentMaterialCutsCoupleIndex();
    const G4PhysicsVector* aVector = nullptr;
    const G4PhysicsVector* bVector = nullptr;
    if(midx < fAsymmetryTable->size())
    {
      aVector = (*fAsymmetryTable)(midx);
    }
    if(midx < fTransverseAsymmetryTable->size())
    {
      bVector = (*fTransverseAsymmetryTable)(midx);
    }
    if(aVector && bVector)
    {
      G4double lAsymmetry = aVector->Value(energy);
      G4double tAsymmetry = bVector->Value(energy);
      G4double polZZ      = polarization.z() * (volPolarization * direction0);
      G4double polXX =
        polarization.x() *
        (volPolarization * G4PolarizationHelper::GetParticleFrameX(direction0));
      G4double polYY =
        polarization.y() *
        (volPolarization * G4PolarizationHelper::GetParticleFrameY(direction0));

      factor /= (1. + polZZ * lAsymmetry + (polXX + polYY) * tAsymmetry);

      if(verboseLevel >= 2)
      {
        G4cout << " Asymmetry:     " << lAsymmetry << ", " << tAsymmetry
               << G4endl;
        G4cout << " PolProduct:    " << polXX << ", " << polYY << ", " << polZZ
               << G4endl;
        G4cout << " Factor:        " << factor << G4endl;
      }
    }
    else
    {
      G4ExceptionDescription ed;
      ed << "Problem with asymmetry tables: material index " << midx
         << " is out of range or tables are not filled";
      G4Exception("G4PolarizedIonisation::ComputeSaturationFactor", "em0048",
                  JustWarning, ed, "");
    }
  }
  return factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedIonisation::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  // *** build DEDX and (unpolarized) cross section tables
  G4VEnergyLossProcess::BuildPhysicsTable(part);
  G4bool master = true;
  const G4PolarizedIonisation* masterProcess =
    static_cast<const G4PolarizedIonisation*>(GetMasterProcess());
  if(masterProcess && masterProcess != this)
  {
    master = false;
  }
  if(master)
  {
    BuildAsymmetryTables(part);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedIonisation::BuildAsymmetryTables(
  const G4ParticleDefinition& part)
{
  // cleanup old, initialise new table
  CleanTables();
  fAsymmetryTable = G4PhysicsTableHelper::PreparePhysicsTable(fAsymmetryTable);
  fTransverseAsymmetryTable =
    G4PhysicsTableHelper::PreparePhysicsTable(fTransverseAsymmetryTable);

  const G4ProductionCutsTable* theCoupleTable =
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

  for(G4int j = 0; j < numOfCouples; ++j)
  {
    // get cut value
    const G4MaterialCutsCouple* couple =
      theCoupleTable->GetMaterialCutsCouple(j);

    G4double cut = (*theCoupleTable->GetEnergyCutsVector(1))[j];

    // create physics vectors then fill it (same parameters as lambda vector)
    G4PhysicsVector* ptrVectorA = LambdaPhysicsVector(couple, cut);
    G4PhysicsVector* ptrVectorB = LambdaPhysicsVector(couple, cut);
    std::size_t bins            = ptrVectorA->GetVectorLength();

    for(std::size_t i = 0; i < bins; ++i)
    {
      G4double lowEdgeEnergy = ptrVectorA->Energy(i);
      G4double tasm          = 0.;
      G4double asym = ComputeAsymmetry(lowEdgeEnergy, couple, part, cut, tasm);
      ptrVectorA->PutValue(i, asym);
      ptrVectorB->PutValue(i, tasm);
    }
    fAsymmetryTable->insertAt(j, ptrVectorA);
    fTransverseAsymmetryTable->insertAt(j, ptrVectorB);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4PolarizedIonisation::ComputeAsymmetry(
  G4double energy, const G4MaterialCutsCouple* couple,
  const G4ParticleDefinition& aParticle, G4double cut, G4double& tAsymmetry)
{
  G4double lAsymmetry = 0.0;
  tAsymmetry          = 0.0;
  if(fIsElectron)
  {
    lAsymmetry = tAsymmetry = -1.0;
  }

  // calculate polarized cross section
  G4ThreeVector targetPolarization = G4ThreeVector(0., 0., 1.);
  fEmModel->SetTargetPolarization(targetPolarization);
  fEmModel->SetBeamPolarization(targetPolarization);
  G4double sigma2 =
    fEmModel->CrossSection(couple, &aParticle, energy, cut, energy);

  // calculate transversely polarized cross section
  targetPolarization = G4ThreeVector(1., 0., 0.);
  fEmModel->SetTargetPolarization(targetPolarization);
  fEmModel->SetBeamPolarization(targetPolarization);
  G4double sigma3 =
    fEmModel->CrossSection(couple, &aParticle, energy, cut, energy);

  // calculate unpolarized cross section
  targetPolarization = G4ThreeVector();
  fEmModel->SetTargetPolarization(targetPolarization);
  fEmModel->SetBeamPolarization(targetPolarization);
  G4double sigma0 =
    fEmModel->CrossSection(couple, &aParticle, energy, cut, energy);
  // determine asymmetries
  if(sigma0 > 0.)
  {
    lAsymmetry = sigma2 / sigma0 - 1.;
    tAsymmetry = sigma3 / sigma0 - 1.;
  }
  if(std::fabs(lAsymmetry) > 1.)
  {
    G4ExceptionDescription ed;
    ed << "G4PolarizedIonisation::ComputeAsymmetry : E(MeV)= " << energy
       << " lAsymmetry= " << lAsymmetry << " (" << std::fabs(lAsymmetry) - 1.
       << ")";
    G4Exception("G4PolarizedIonisation::ComputeAsymmetry", "pol002",
                JustWarning, ed);
  }
  if(std::fabs(tAsymmetry) > 1.)
  {
    G4ExceptionDescription ed;
    ed << "G4PolarizedIonisation::ComputeAsymmetry : E(MeV)= " << energy
       << " tAsymmetry= " << tAsymmetry << " (" << std::fabs(tAsymmetry) - 1.
       << ")";
    G4Exception("G4PolarizedIonisation::ComputeAsymmetry", "pol003",
                JustWarning, ed);
  }
  return lAsymmetry;
}
