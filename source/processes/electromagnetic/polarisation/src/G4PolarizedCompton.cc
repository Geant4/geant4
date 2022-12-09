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
// File name:     G4PolarizedCompton
//
// Author:        Andreas Schaelicke
//                based on code by Michel Maire / Vladimir IVANTCHENKO
//
// Class description
//   modified version respecting media and beam polarization
//   using the stokes formalism

#include "G4PolarizedCompton.hh"

#include "G4Electron.hh"
#include "G4EmParameters.hh"
#include "G4KleinNishinaCompton.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4PolarizationManager.hh"
#include "G4PolarizedComptonModel.hh"
#include "G4ProductionCutsTable.hh"
#include "G4StokesVector.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PhysicsTable* G4PolarizedCompton::theAsymmetryTable = nullptr;

G4PolarizedCompton::G4PolarizedCompton(const G4String& processName,
                                       G4ProcessType type)
  : G4VEmProcess(processName, type)
  , fType(10)
  , fBuildAsymmetryTable(true)
  , fUseAsymmetryTable(true)
  , fIsInitialised(false)
{
  SetStartFromNullFlag(true);
  SetBuildTableFlag(true);
  SetSecondaryParticle(G4Electron::Electron());
  SetProcessSubType(fComptonScattering);
  SetMinKinEnergyPrim(1. * MeV);
  SetSplineFlag(true);
  fEmModel = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PolarizedCompton::~G4PolarizedCompton() { CleanTable(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedCompton::ProcessDescription(std::ostream& out) const
{
  out << "Polarized model for Compton scattering.\n";

  G4VEmProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedCompton::CleanTable()
{
  if(theAsymmetryTable)
  {
    theAsymmetryTable->clearAndDestroy();
    delete theAsymmetryTable;
    theAsymmetryTable = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4bool G4PolarizedCompton::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PolarizedCompton::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!fIsInitialised)
  {
    fIsInitialised = true;
    if(0 == fType)
    {
      if(nullptr == EmModel(0))
      {
        SetEmModel(new G4KleinNishinaCompton());
      }
    }
    else
    {
      fEmModel = new G4PolarizedComptonModel();
      SetEmModel(fEmModel);
    }
    G4EmParameters* param = G4EmParameters::Instance();
    EmModel(0)->SetLowEnergyLimit(param->MinKinEnergy());
    EmModel(0)->SetHighEnergyLimit(param->MaxKinEnergy());
    AddEmModel(1, EmModel(0));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedCompton::SetModel(const G4String& ss)
{
  if(ss == "Klein-Nishina")
  {
    fType = 0;
  }
  if(ss == "Polarized-Compton")
  {
    fType = 10;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4PolarizedCompton::GetMeanFreePath(const G4Track& aTrack,
                                             G4double previousStepSize,
                                             G4ForceCondition* condition)
{
  // *** get unploarised mean free path from lambda table ***
  G4double mfp =
    G4VEmProcess::GetMeanFreePath(aTrack, previousStepSize, condition);

  if(theAsymmetryTable && fUseAsymmetryTable && mfp < DBL_MAX)
  {
    mfp *= ComputeSaturationFactor(aTrack);
  }
  if(verboseLevel >= 2)
  {
    G4cout << "G4PolarizedCompton::MeanFreePath:  " << mfp / mm << " mm "
           << G4endl;
  }
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4PolarizedCompton::PostStepGetPhysicalInteractionLength(
  const G4Track& aTrack, G4double previousStepSize, G4ForceCondition* condition)
{
  // save previous values
  G4double nLength = theNumberOfInteractionLengthLeft;
  G4double iLength = currentInteractionLength;

  // *** compute unpolarized step limit ***
  // this changes theNumberOfInteractionLengthLeft and currentInteractionLength
  G4double x = G4VEmProcess::PostStepGetPhysicalInteractionLength(
    aTrack, previousStepSize, condition);
  G4double x0      = x;
  G4double satFact = 1.0;

  // *** add corrections on polarisation ***
  if(theAsymmetryTable && fUseAsymmetryTable && x < DBL_MAX)
  {
    satFact            = ComputeSaturationFactor(aTrack);
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
    G4cout << "G4PolarizedCompton::PostStepGPIL: " << std::setprecision(8)
           << x / mm << " mm;" << G4endl
           << "               unpolarized value: " << std::setprecision(8)
           << x0 / mm << " mm." << G4endl;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4PolarizedCompton::ComputeSaturationFactor(const G4Track& aTrack)
{
  G4double factor = 1.0;

  // *** get asymmetry, if target is polarized ***
  const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
  const G4double GammaEnergy             = aDynamicGamma->GetKineticEnergy();
  const G4StokesVector GammaPolarization =
    G4StokesVector(aTrack.GetPolarization());
  const G4ParticleMomentum GammaDirection0 =
    aDynamicGamma->GetMomentumDirection();

  const G4Material* aMaterial = aTrack.GetMaterial();
  G4VPhysicalVolume* aPVolume = aTrack.GetVolume();
  G4LogicalVolume* aLVolume   = aPVolume->GetLogicalVolume();

  G4PolarizationManager* polarizationManager =
    G4PolarizationManager::GetInstance();

  const G4bool VolumeIsPolarized = polarizationManager->IsPolarized(aLVolume);
  G4StokesVector ElectronPolarization =
    polarizationManager->GetVolumePolarization(aLVolume);

  if(VolumeIsPolarized)
  {
    if(verboseLevel >= 2)
    {
      G4cout << "G4PolarizedCompton::ComputeSaturationFactor: " << G4endl;
      G4cout << " Mom " << GammaDirection0 << G4endl;
      G4cout << " Polarization " << GammaPolarization << G4endl;
      G4cout << " MaterialPol. " << ElectronPolarization << G4endl;
      G4cout << " Phys. Volume " << aPVolume->GetName() << G4endl;
      G4cout << " Log. Volume  " << aLVolume->GetName() << G4endl;
      G4cout << " Material     " << aMaterial << G4endl;
    }

    std::size_t midx               = CurrentMaterialCutsCoupleIndex();
    const G4PhysicsVector* aVector = nullptr;
    if(midx < theAsymmetryTable->size())
    {
      aVector = (*theAsymmetryTable)(midx);
    }
    if(aVector)
    {
      G4double asymmetry = aVector->Value(GammaEnergy);

      //  we have to determine angle between particle motion
      //  and target polarisation here
      //      circ pol * Vec(ElectronPol)*Vec(PhotonMomentum)
      //  both vectors in global reference frame

      G4double pol        = ElectronPolarization * GammaDirection0;
      G4double polProduct = GammaPolarization.p3() * pol;
      factor /= (1. + polProduct * asymmetry);
      if(verboseLevel >= 2)
      {
        G4cout << " Asymmetry:     " << asymmetry << G4endl;
        G4cout << " PolProduct:    " << polProduct << G4endl;
        G4cout << " Factor:        " << factor << G4endl;
      }
    }
    else
    {
      G4ExceptionDescription ed;
      ed << "Problem with asymmetry table: material index " << midx
         << " is out of range or the table is not filled";
      G4Exception("G4PolarizedComptonModel::ComputeSaturationFactor", "em0048",
                  JustWarning, ed, "");
    }
  }
  return factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedCompton::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  // *** build (unpolarized) cross section tables (Lambda)
  G4VEmProcess::BuildPhysicsTable(part);
  if(fBuildAsymmetryTable && fEmModel)
  {
    G4bool isMaster = true;
    const G4PolarizedCompton* masterProcess =
      static_cast<const G4PolarizedCompton*>(GetMasterProcess());
    if(masterProcess && masterProcess != this)
    {
      isMaster = false;
    }
    if(isMaster)
    {
      BuildAsymmetryTable(part);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedCompton::BuildAsymmetryTable(const G4ParticleDefinition& part)
{
  // cleanup old, initialise new table
  CleanTable();
  theAsymmetryTable =
    G4PhysicsTableHelper::PreparePhysicsTable(theAsymmetryTable);

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable =
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();
  if(!theAsymmetryTable)
  {
    return;
  }
  G4int nbins                 = LambdaBinning();
  G4double emin               = MinKinEnergy();
  G4double emax               = MaxKinEnergy();
  G4PhysicsLogVector* aVector = nullptr;
  G4PhysicsLogVector* bVector = nullptr;

  for(G4int i = 0; i < numOfCouples; ++i)
  {
    if(theAsymmetryTable->GetFlag(i))
    {
      // create physics vector and fill it
      const G4MaterialCutsCouple* couple =
        theCoupleTable->GetMaterialCutsCouple(i);
      // use same parameters as for lambda
      if(!aVector)
      {
        aVector = new G4PhysicsLogVector(emin, emax, nbins, true);
        bVector = aVector;
      }
      else
      {
        bVector = new G4PhysicsLogVector(*aVector);
      }

      for(G4int j = 0; j <= nbins; ++j)
      {
        G4double energy = bVector->Energy(j);
        G4double tasm   = 0.;
        G4double asym   = ComputeAsymmetry(energy, couple, part, 0., tasm);
        bVector->PutValue(j, asym);
      }
      bVector->FillSecondDerivatives();
      G4PhysicsTableHelper::SetPhysicsVector(theAsymmetryTable, i, bVector);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4PolarizedCompton::ComputeAsymmetry(
  G4double energy, const G4MaterialCutsCouple* couple,
  const G4ParticleDefinition& aParticle, G4double cut, G4double& tAsymmetry)
{
  G4double lAsymmetry = 0.0;
  tAsymmetry          = 0;

  // calculate polarized cross section
  G4ThreeVector thePolarization = G4ThreeVector(0., 0., 1.);
  fEmModel->SetTargetPolarization(thePolarization);
  fEmModel->SetBeamPolarization(thePolarization);
  G4double sigma2 =
    fEmModel->CrossSection(couple, &aParticle, energy, cut, energy);

  // calculate unpolarized cross section
  thePolarization = G4ThreeVector();
  fEmModel->SetTargetPolarization(thePolarization);
  fEmModel->SetBeamPolarization(thePolarization);
  G4double sigma0 =
    fEmModel->CrossSection(couple, &aParticle, energy, cut, energy);

  // determine asymmetries
  if(sigma0 > 0.)
  {
    lAsymmetry = sigma2 / sigma0 - 1.;
  }
  return lAsymmetry;
}
