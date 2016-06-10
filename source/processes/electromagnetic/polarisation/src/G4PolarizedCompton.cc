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
// $Id: G4PolarizedCompton.cc 85018 2014-10-23 09:51:37Z gcosmo $
// 
//
// File name:     G4PolarizedCompton
//
// Author:        Andreas Schaelicke
//                based on code by Michel Maire / Vladimir IVANTCHENKO
// Class description
//
// modified version respecting media and beam polarization
//     using the stokes formalism
//
// Creation date: 01.05.2005
//
// Modifications:
//
// 01-01-05, include polarization description (A.Stahl)
// 01-01-05, create asymmetry table and determine interactionlength (A.Stahl)
// 01-05-05, update handling of media polarization (A.Schalicke)
// 01-05-05, update polarized differential cross section (A.Schalicke)
// 20-05-05, added polarization transfer (A.Schalicke)
// 10-06-05, transformation between different reference frames (A.Schalicke)
// 17-10-05, correct reference frame dependence in GetMeanFreePath (A.Schalicke)
// 26-07-06, cross section recalculated (P.Starovoitov)
// 09-08-06, make it work under current geant4 release (A.Schalicke)
// 11-06-07, add PostStepGetPhysicalInteractionLength (A.Schalicke)
// -----------------------------------------------------------------------------


#include "G4PolarizedCompton.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"

#include "G4StokesVector.hh"
#include "G4PolarizationManager.hh"
#include "G4PolarizedComptonModel.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4KleinNishinaCompton.hh"
#include "G4PolarizedComptonModel.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PhysicsTable* G4PolarizedCompton::theAsymmetryTable = 0;

G4PolarizedCompton::G4PolarizedCompton(const G4String& processName,
  G4ProcessType type):
  G4VEmProcess (processName, type),
  buildAsymmetryTable(true),
  useAsymmetryTable(true),
  isInitialised(false),
  mType(10)
{
  SetStartFromNullFlag(true);
  SetBuildTableFlag(true);
  SetSecondaryParticle(G4Electron::Electron());
  SetProcessSubType(fComptonScattering);
  SetMinKinEnergyPrim(1*MeV);
  SetSplineFlag(true);
  emModel = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4PolarizedCompton::~G4PolarizedCompton()
{
  delete theAsymmetryTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PolarizedCompton::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!isInitialised) {
    isInitialised = true;
    if(0 == mType) {
      if(!EmModel(1)) { SetEmModel(new G4KleinNishinaCompton(), 1); }
    } else {
      emModel = new G4PolarizedComptonModel();
      SetEmModel(emModel, 1); 
    }
    G4EmParameters* param = G4EmParameters::Instance();
    EmModel(1)->SetLowEnergyLimit(param->MinKinEnergy());
    EmModel(1)->SetHighEnergyLimit(param->MaxKinEnergy());
    AddEmModel(1, EmModel(1));
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PolarizedCompton::PrintInfo()
{
  G4cout << " Total cross sections has a good parametrisation"
         << " from 10 KeV to (100/Z) GeV" 
         << "\n      Sampling according " <<  EmModel(1)->GetName() << " model" 
	 << G4endl;
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PolarizedCompton::SetModel(const G4String& ss)
{
  if(ss == "Klein-Nishina")     { mType = 0; }
  if(ss == "Polarized-Compton") { mType = 10; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4PolarizedCompton::GetMeanFreePath(
				   const G4Track& aTrack,
				   G4double   previousStepSize,
				   G4ForceCondition* condition)
{
  // *** get unploarised mean free path from lambda table ***
  G4double mfp = G4VEmProcess::GetMeanFreePath(aTrack, previousStepSize, condition);

  if (theAsymmetryTable && useAsymmetryTable) {
    // *** get asymmetry, if target is polarized ***
    const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
    const G4double GammaEnergy = aDynamicGamma->GetKineticEnergy();
    const G4StokesVector GammaPolarization = aTrack.GetPolarization();
    const G4ParticleMomentum GammaDirection0 = aDynamicGamma->GetMomentumDirection();

    G4Material*         aMaterial = aTrack.GetMaterial();
    G4VPhysicalVolume*  aPVolume  = aTrack.GetVolume();
    G4LogicalVolume*    aLVolume  = aPVolume->GetLogicalVolume();

    //   G4Material* bMaterial = aLVolume->GetMaterial();
    G4PolarizationManager * polarizationManger = G4PolarizationManager::GetInstance();

    G4bool VolumeIsPolarized = polarizationManger->IsPolarized(aLVolume);
    G4StokesVector ElectronPolarization = polarizationManger->GetVolumePolarization(aLVolume);
     
    if (!VolumeIsPolarized || mfp == DBL_MAX) return mfp;
     
    if (verboseLevel>=2) {

      G4cout << " Mom " << GammaDirection0  << G4endl;
      G4cout << " Polarization " << GammaPolarization  << G4endl;
      G4cout << " MaterialPol. " << ElectronPolarization  << G4endl;
      G4cout << " Phys. Volume " << aPVolume->GetName() << G4endl;
      G4cout << " Log. Volume  " << aLVolume->GetName() << G4endl;
      G4cout << " Material     " << aMaterial          << G4endl;
    }

    G4int midx= CurrentMaterialCutsCoupleIndex();
    G4PhysicsVector * aVector=(*theAsymmetryTable)(midx);
     
    G4double asymmetry=0;
    if (aVector) {
      asymmetry = aVector->Value(GammaEnergy);
    } else {
      G4cout << " MaterialIndex     " << midx << " is out of range \n";
      asymmetry=0;
    }

    //  we have to determine angle between particle motion 
    //  and target polarisation here  
    //      circ pol * Vec(ElectronPol)*Vec(PhotonMomentum)
    //  both vectors in global reference frame
     
    G4double pol=ElectronPolarization*GammaDirection0;
     
    G4double polProduct = GammaPolarization.p3() * pol;
    mfp *= 1. / ( 1. + polProduct * asymmetry );

    if (verboseLevel>=2) {
      G4cout << " MeanFreePath:  " << mfp / mm << " mm " << G4endl;
      G4cout << " Asymmetry:     " << asymmetry          << G4endl;
      G4cout << " PolProduct:    " << polProduct         << G4endl;
    }
  }

  return mfp;
}

G4double G4PolarizedCompton::PostStepGetPhysicalInteractionLength(
				   const G4Track& aTrack,
				   G4double   previousStepSize,
				   G4ForceCondition* condition)
{
  // *** get unploarised mean free path from lambda table ***
  G4double mfp = G4VEmProcess::PostStepGetPhysicalInteractionLength(aTrack, previousStepSize, condition);


   if (theAsymmetryTable && useAsymmetryTable) {
     // *** get asymmetry, if target is polarized ***
     const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
     const G4double GammaEnergy = aDynamicGamma->GetKineticEnergy();
     const G4StokesVector GammaPolarization = aTrack.GetPolarization();
     const G4ParticleMomentum GammaDirection0 = aDynamicGamma->GetMomentumDirection();

     G4Material*         aMaterial = aTrack.GetMaterial();
     G4VPhysicalVolume*  aPVolume  = aTrack.GetVolume();
     G4LogicalVolume*    aLVolume  = aPVolume->GetLogicalVolume();

     //   G4Material* bMaterial = aLVolume->GetMaterial();
     G4PolarizationManager * polarizationManger = G4PolarizationManager::GetInstance();

     const G4bool VolumeIsPolarized = polarizationManger->IsPolarized(aLVolume);
     G4StokesVector ElectronPolarization = polarizationManger->GetVolumePolarization(aLVolume);

     if (!VolumeIsPolarized || mfp == DBL_MAX) return mfp;
     
     if (verboseLevel>=2) {

       G4cout << " Mom " << GammaDirection0  << G4endl;
       G4cout << " Polarization " << GammaPolarization  << G4endl;
       G4cout << " MaterialPol. " << ElectronPolarization  << G4endl;
       G4cout << " Phys. Volume " << aPVolume->GetName() << G4endl;
       G4cout << " Log. Volume  " << aLVolume->GetName() << G4endl;
       G4cout << " Material     " << aMaterial          << G4endl;
     }

     G4int midx= CurrentMaterialCutsCoupleIndex();
     G4PhysicsVector * aVector=(*theAsymmetryTable)(midx);
     
     G4double asymmetry=0;
     if (aVector) {
       asymmetry = aVector->Value(GammaEnergy);
     } else {
       G4cout << " MaterialIndex     " << midx << " is out of range \n";
       asymmetry=0;
     }

     //  we have to determine angle between particle motion 
     //  and target polarisation here  
     //      circ pol * Vec(ElectronPol)*Vec(PhotonMomentum)
     //  both vectors in global reference frame
     
     G4double pol=ElectronPolarization*GammaDirection0;
     
     G4double polProduct = GammaPolarization.p3() * pol;
     mfp *= 1. / ( 1. + polProduct * asymmetry );

     if (verboseLevel>=2) {
       G4cout << " MeanFreePath:  " << mfp / mm << " mm " << G4endl;
       G4cout << " Asymmetry:     " << asymmetry          << G4endl;
       G4cout << " PolProduct:    " << polProduct         << G4endl;
     }
   }

   return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4PolarizedCompton::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  G4VEmProcess::PreparePhysicsTable(part);

  if(buildAsymmetryTable && emModel) {
    G4bool isMaster = true;
    const G4PolarizedCompton* masterProcess = 
      static_cast<const G4PolarizedCompton*>(GetMasterProcess());
    if(masterProcess && masterProcess != this) { isMaster = false; }
    if(isMaster) {
      theAsymmetryTable = 
	G4PhysicsTableHelper::PreparePhysicsTable(theAsymmetryTable);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void G4PolarizedCompton::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  // *** build (unpolarized) cross section tables (Lambda)
  G4VEmProcess::BuildPhysicsTable(part);
  if(buildAsymmetryTable && emModel) { 
    G4bool isMaster = true;
    const G4PolarizedCompton* masterProcess = 
      static_cast<const G4PolarizedCompton*>(GetMasterProcess());
    if(masterProcess && masterProcess != this) { isMaster = false; }
    if(isMaster) { BuildAsymmetryTable(part); }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void G4PolarizedCompton::BuildAsymmetryTable(const G4ParticleDefinition& part)
{
  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  if(!theAsymmetryTable) { return; }
  G4int nbins = LambdaBinning();
  G4double emin = MinKinEnergy();
  G4double emax = MaxKinEnergy();
  G4PhysicsLogVector* aVector = 0;
  G4PhysicsLogVector* bVector = 0;

  for(size_t i=0; i<numOfCouples; ++i) {
    if (theAsymmetryTable->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
      // use same parameters as for lambda
      if(!aVector) { 
	aVector = new G4PhysicsLogVector(emin, emax, nbins); 
        aVector->SetSpline(true);
        bVector = aVector;
      } else {
	bVector = new G4PhysicsLogVector(*aVector);
      }

      for (G4int j = 0; j <= nbins; ++j ) {
	G4double energy = bVector->Energy(j);
	G4double tasm=0.;
	G4double asym = ComputeAsymmetry(energy, couple, part, 0., tasm);
	bVector->PutValue(j,asym);
      }

      G4PhysicsTableHelper::SetPhysicsVector(theAsymmetryTable, i, bVector);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double G4PolarizedCompton::ComputeAsymmetry(G4double energy,
					      const G4MaterialCutsCouple* couple,
					      const G4ParticleDefinition& aParticle,
					      G4double cut,
					      G4double & tAsymmetry)
{
  G4double lAsymmetry = 0.0;
  tAsymmetry=0;

  //
  // calculate polarized cross section
  //
  G4ThreeVector thePolarization=G4ThreeVector(0.,0.,1.);
  emModel->SetTargetPolarization(thePolarization);
  emModel->SetBeamPolarization(thePolarization);
  G4double sigma2=emModel->CrossSection(couple,&aParticle,energy,cut,energy);

  //
  // calculate unpolarized cross section
  //
  thePolarization=G4ThreeVector();
  emModel->SetTargetPolarization(thePolarization);
  emModel->SetBeamPolarization(thePolarization);
  G4double sigma0=emModel->CrossSection(couple,&aParticle,energy,cut,energy);

  // determine assymmetries
  if (sigma0 > 0.) {
    lAsymmetry = sigma2/sigma0-1.;
  }
  return lAsymmetry;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
