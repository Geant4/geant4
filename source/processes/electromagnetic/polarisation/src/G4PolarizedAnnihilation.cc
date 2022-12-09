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
// File name:     G4PolarizedAnnihilation
//
// Author:  A. Schaelicke on base of Vladimir Ivanchenko / Michel Maire code
//
// Class Description:
//   Polarized process of e+ annihilation into 2 gammas

#include "G4PolarizedAnnihilation.hh"

#include "G4DynamicParticle.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4PhysicsVector.hh"
#include "G4PolarizationHelper.hh"
#include "G4PolarizationManager.hh"
#include "G4PolarizedAnnihilationModel.hh"
#include "G4ProductionCutsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4StokesVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4PolarizedAnnihilation::G4PolarizedAnnihilation(const G4String& name)
  : G4eplusAnnihilation(name)
  , fAsymmetryTable(nullptr)
  , fTransverseAsymmetryTable(nullptr)
{
  fEmModel = new G4PolarizedAnnihilationModel();
  SetEmModel(fEmModel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4PolarizedAnnihilation::~G4PolarizedAnnihilation() { CleanTables(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedAnnihilation::CleanTables()
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
G4double G4PolarizedAnnihilation::GetMeanFreePath(const G4Track& track,
                                                  G4double previousStepSize,
                                                  G4ForceCondition* condition)
{
  G4double mfp =
    G4VEmProcess::GetMeanFreePath(track, previousStepSize, condition);

  if(nullptr != fAsymmetryTable && nullptr != fTransverseAsymmetryTable && mfp < DBL_MAX)
  {
    mfp *= ComputeSaturationFactor(track);
  }
  if(verboseLevel >= 2)
  {
    G4cout << "G4PolarizedAnnihilation::MeanFreePath:  " << mfp / mm << " mm "
           << G4endl;
  }
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4PolarizedAnnihilation::PostStepGetPhysicalInteractionLength(
  const G4Track& track, G4double previousStepSize, G4ForceCondition* condition)
{
  // save previous values
  G4double nLength = theNumberOfInteractionLengthLeft;
  G4double iLength = currentInteractionLength;

  // *** compute unpolarized step limit ***
  // this changes theNumberOfInteractionLengthLeft and currentInteractionLength
  G4double x = G4VEmProcess::PostStepGetPhysicalInteractionLength(
    track, previousStepSize, condition);
  G4double x0      = x;
  G4double satFact = 1.0;

  // *** add corrections on polarisation ***
  if(nullptr != fAsymmetryTable && nullptr != fTransverseAsymmetryTable && x < DBL_MAX)
  {
    satFact            = ComputeSaturationFactor(track);
    G4double curLength = currentInteractionLength * satFact;
    G4double prvLength = iLength * satFact;
    if(nLength > 0.0)
    {
      theNumberOfInteractionLengthLeft =
        std::max(nLength - previousStepSize / prvLength, 0.0);
    }
    x = theNumberOfInteractionLengthLeft * curLength;
  }
  if(verboseLevel >= 2)
  {
    G4cout << "G4PolarizedAnnihilation::PostStepGPIL: " << std::setprecision(8)
           << x / mm << " mm;" << G4endl
           << "                         unpolarized value: "
           << std::setprecision(8) << x0 / mm << " mm." << G4endl;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4PolarizedAnnihilation::ComputeSaturationFactor(const G4Track& track)
{
  const G4Material* aMaterial = track.GetMaterial();
  G4VPhysicalVolume* aPVolume = track.GetVolume();
  G4LogicalVolume* aLVolume   = aPVolume->GetLogicalVolume();

  G4PolarizationManager* polarizationManager =
    G4PolarizationManager::GetInstance();

  const G4bool volumeIsPolarized = polarizationManager->IsPolarized(aLVolume);
  G4StokesVector electronPolarization =
    polarizationManager->GetVolumePolarization(aLVolume);

  G4double factor = 1.0;

  if(volumeIsPolarized)
  {
    // *** get asymmetry, if target is polarized ***
    const G4DynamicParticle* aDynamicPositron = track.GetDynamicParticle();
    const G4double positronEnergy = aDynamicPositron->GetKineticEnergy();
    const G4StokesVector positronPolarization =
      G4StokesVector(track.GetPolarization());
    const G4ParticleMomentum positronDirection0 =
      aDynamicPositron->GetMomentumDirection();

    if(verboseLevel >= 2)
    {
      G4cout << "G4PolarizedAnnihilation::ComputeSaturationFactor: " << G4endl;
      G4cout << " Mom " << positronDirection0 << G4endl;
      G4cout << " Polarization " << positronPolarization << G4endl;
      G4cout << " MaterialPol. " << electronPolarization << G4endl;
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
      G4double lAsymmetry = aVector->Value(positronEnergy);
      G4double tAsymmetry = bVector->Value(positronEnergy);
      G4double polZZ =
        positronPolarization.z() * (electronPolarization * positronDirection0);
      G4double polXX =
        positronPolarization.x() *
        (electronPolarization *
         G4PolarizationHelper::GetParticleFrameX(positronDirection0));
      G4double polYY =
        positronPolarization.y() *
        (electronPolarization *
         G4PolarizationHelper::GetParticleFrameY(positronDirection0));

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
      G4Exception("G4PolarizedAnnihilation::ComputeSaturationFactor", "em0048",
                  JustWarning, ed, "");
    }
  }
  return factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedAnnihilation::BuildPhysicsTable(
  const G4ParticleDefinition& part)
{
  G4VEmProcess::BuildPhysicsTable(part);
  if(isTheMaster)
  {
    BuildAsymmetryTables(part);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedAnnihilation::BuildAsymmetryTables(
  const G4ParticleDefinition& part)
{
  // cleanup old, initialise new table
  CleanTables();
  fAsymmetryTable = G4PhysicsTableHelper::PreparePhysicsTable(fAsymmetryTable);
  fTransverseAsymmetryTable =
    G4PhysicsTableHelper::PreparePhysicsTable(fTransverseAsymmetryTable);
  if(nullptr == fAsymmetryTable) return;

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable =
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();
  for(G4int i = 0; i < numOfCouples; ++i)
  {
    if(fAsymmetryTable->GetFlag(i))
    {
      // create physics vector and fill it
      const G4MaterialCutsCouple* couple =
        theCoupleTable->GetMaterialCutsCouple(i);

      // use same parameters as for lambda
      G4PhysicsVector* aVector = LambdaPhysicsVector(couple);
      G4PhysicsVector* tVector = LambdaPhysicsVector(couple);
      G4int nn = (G4int)aVector->GetVectorLength();
      for(G4int j = 0; j < nn; ++j)
      {
        G4double energy = aVector->Energy(j);
        G4double tasm   = 0.;
        G4double asym = ComputeAsymmetry(energy, couple, part, 0., tasm);
        aVector->PutValue(j, asym);
        tVector->PutValue(j, tasm);
      }
      if(aVector->GetSpline()) {
	aVector->FillSecondDerivatives();
	tVector->FillSecondDerivatives();
      }
      G4PhysicsTableHelper::SetPhysicsVector(fAsymmetryTable, i, aVector);
      G4PhysicsTableHelper::SetPhysicsVector(fTransverseAsymmetryTable, i,
                                             tVector);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4PolarizedAnnihilation::ComputeAsymmetry(
  G4double energy, const G4MaterialCutsCouple* couple,
  const G4ParticleDefinition& aParticle, G4double cut, G4double& tAsymmetry)
{
  G4double lAsymmetry = 0.0;
  tAsymmetry          = 0.0;

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
  return lAsymmetry;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedAnnihilation::ProcessDescription(std::ostream& out) const
{
  out << "Polarized model for positron annihilation into 2 photons.\n";
}
