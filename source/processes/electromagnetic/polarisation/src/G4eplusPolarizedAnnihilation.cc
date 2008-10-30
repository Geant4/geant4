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
// $Id: G4eplusPolarizedAnnihilation.cc,v 1.8 2008-10-30 22:34:23 schaelic Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eplusPolarizedAnnihilation
//
// Author:        A. Schaelicke on base of Vladimir Ivanchenko / Michel Maire code
//
// Creation date: 02.07.2006
//
// Modifications:
// 26-07-06 modified cross section  (P. Starovoitov)
// 21-08-06 interface updated   (A. Schaelicke)
// 11-06-07, add PostStepGetPhysicalInteractionLength (A.Schalicke)
// 02-10-07, enable AtRest (V.Ivanchenko)
//
//
// Class Description:
//
// Polarized process of e+ annihilation into 2 gammas
//

//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eplusPolarizedAnnihilation.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Gamma.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLogVector.hh"


#include "G4PolarizedAnnihilationModel.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PolarizationManager.hh"
#include "G4PolarizationHelper.hh"
#include "G4StokesVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eplusPolarizedAnnihilation::G4eplusPolarizedAnnihilation(const G4String& name)
  : G4VEmProcess(name), isInitialised(false),
    theAsymmetryTable(NULL),
    theTransverseAsymmetryTable(NULL)
{
  enableAtRestDoIt = true;
  SetProcessSubType(fAnnihilation);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eplusPolarizedAnnihilation::~G4eplusPolarizedAnnihilation()
{
  if (theAsymmetryTable) {
    theAsymmetryTable->clearAndDestroy();
    delete theAsymmetryTable;
  }
  if (theTransverseAsymmetryTable) {
    theTransverseAsymmetryTable->clearAndDestroy();
    delete theTransverseAsymmetryTable;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusPolarizedAnnihilation::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!isInitialised) {
    isInitialised = true;
    //    SetVerboseLevel(3);
    SetBuildTableFlag(true);
    SetStartFromNullFlag(false);
    SetSecondaryParticle(G4Gamma::Gamma());
    G4double emin = 0.1*keV;
    G4double emax = 100.*TeV;
    SetLambdaBinning(120);
    SetMinKinEnergy(emin);
    SetMaxKinEnergy(emax);
    emModel = new G4PolarizedAnnihilationModel();
    emModel->SetLowEnergyLimit(emin);
    emModel->SetHighEnergyLimit(emax);
    AddEmModel(1, emModel);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  // for polarization

G4double G4eplusPolarizedAnnihilation::GetMeanFreePath(const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition)
{
  G4double mfp = G4VEmProcess::GetMeanFreePath(track, previousStepSize, condition);

  if (theAsymmetryTable) {

    G4Material*         aMaterial = track.GetMaterial();
    G4VPhysicalVolume*  aPVolume  = track.GetVolume();
    G4LogicalVolume*    aLVolume  = aPVolume->GetLogicalVolume();
    
    //   G4Material* bMaterial = aLVolume->GetMaterial();
    G4PolarizationManager * polarizationManger = G4PolarizationManager::GetInstance();
    
    const G4bool volumeIsPolarized = polarizationManger->IsPolarized(aLVolume);
    G4StokesVector electronPolarization = polarizationManger->GetVolumePolarization(aLVolume);

    if (!volumeIsPolarized || mfp == DBL_MAX) return mfp;
     
    // *** get asymmetry, if target is polarized ***
    const G4DynamicParticle* aDynamicPositron = track.GetDynamicParticle();
    const G4double positronEnergy = aDynamicPositron->GetKineticEnergy();
    const G4StokesVector positronPolarization = track.GetPolarization();
    const G4ParticleMomentum positronDirection0 = aDynamicPositron->GetMomentumDirection();

    if (verboseLevel>=2) {
      
      G4cout << " Mom " << positronDirection0  << G4endl;
      G4cout << " Polarization " << positronPolarization  << G4endl;
      G4cout << " MaterialPol. " << electronPolarization  << G4endl;
      G4cout << " Phys. Volume " << aPVolume->GetName() << G4endl;
      G4cout << " Log. Volume  " << aLVolume->GetName() << G4endl;
      G4cout << " Material     " << aMaterial          << G4endl;
    }
    
    G4bool isOutRange;
    G4int idx= CurrentMaterialCutsCoupleIndex();
    G4double lAsymmetry = (*theAsymmetryTable)(idx)->
                                  GetValue(positronEnergy, isOutRange);
    G4double tAsymmetry = (*theTransverseAsymmetryTable)(idx)->
                                  GetValue(positronEnergy, isOutRange);

    G4double polZZ = positronPolarization.z()*
      electronPolarization*positronDirection0;
    G4double polXX = positronPolarization.x()*
      electronPolarization*G4PolarizationHelper::GetParticleFrameX(positronDirection0);
    G4double polYY = positronPolarization.y()*
      electronPolarization*G4PolarizationHelper::GetParticleFrameY(positronDirection0);

    G4double impact = 1. + polZZ*lAsymmetry + (polXX + polYY)*tAsymmetry;

    mfp *= 1. / impact;

    if (verboseLevel>=2) {
      G4cout << " MeanFreePath:  " << mfp / mm << " mm " << G4endl;
      G4cout << " Asymmetry:     " << lAsymmetry << ", " << tAsymmetry  << G4endl;
      G4cout << " PolProduct:    " << polXX << ", " << polYY << ", " << polZZ << G4endl;
    }
  }
  
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eplusPolarizedAnnihilation::PostStepGetPhysicalInteractionLength(
                              const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition)
{
  G4double mfp = G4VEmProcess::PostStepGetPhysicalInteractionLength(track, previousStepSize, condition);

  if (theAsymmetryTable) {

    G4Material*         aMaterial = track.GetMaterial();
    G4VPhysicalVolume*  aPVolume  = track.GetVolume();
    G4LogicalVolume*    aLVolume  = aPVolume->GetLogicalVolume();
    
    //   G4Material* bMaterial = aLVolume->GetMaterial();
    G4PolarizationManager * polarizationManger = G4PolarizationManager::GetInstance();
    
    const G4bool volumeIsPolarized = polarizationManger->IsPolarized(aLVolume);
    G4StokesVector electronPolarization = polarizationManger->GetVolumePolarization(aLVolume);

    if (!volumeIsPolarized || mfp == DBL_MAX) return mfp;
     
    // *** get asymmetry, if target is polarized ***
    const G4DynamicParticle* aDynamicPositron = track.GetDynamicParticle();
    const G4double positronEnergy = aDynamicPositron->GetKineticEnergy();
    const G4StokesVector positronPolarization = track.GetPolarization();
    const G4ParticleMomentum positronDirection0 = aDynamicPositron->GetMomentumDirection();

    if (verboseLevel>=2) {
      
      G4cout << " Mom " << positronDirection0  << G4endl;
      G4cout << " Polarization " << positronPolarization  << G4endl;
      G4cout << " MaterialPol. " << electronPolarization  << G4endl;
      G4cout << " Phys. Volume " << aPVolume->GetName() << G4endl;
      G4cout << " Log. Volume  " << aLVolume->GetName() << G4endl;
      G4cout << " Material     " << aMaterial          << G4endl;
    }
    
    G4bool isOutRange;
    G4int idx= CurrentMaterialCutsCoupleIndex();
    G4double lAsymmetry = (*theAsymmetryTable)(idx)->
                                  GetValue(positronEnergy, isOutRange);
    G4double tAsymmetry = (*theTransverseAsymmetryTable)(idx)->
                                  GetValue(positronEnergy, isOutRange);

    G4double polZZ = positronPolarization.z()*
      electronPolarization*positronDirection0;
    G4double polXX = positronPolarization.x()*
      electronPolarization*G4PolarizationHelper::GetParticleFrameX(positronDirection0);
    G4double polYY = positronPolarization.y()*
      electronPolarization*G4PolarizationHelper::GetParticleFrameY(positronDirection0);

    G4double impact = 1. + polZZ*lAsymmetry + (polXX + polYY)*tAsymmetry;

    mfp *= 1. / impact;

    if (verboseLevel>=2) {
      G4cout << " MeanFreePath:  " << mfp / mm << " mm " << G4endl;
      G4cout << " Asymmetry:     " << lAsymmetry << ", " << tAsymmetry  << G4endl;
      G4cout << " PolProduct:    " << polXX << ", " << polYY << ", " << polZZ << G4endl;
    }
  }
  
  return mfp;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusPolarizedAnnihilation::BuildPhysicsTable(const G4ParticleDefinition& pd) 
{
  G4VEmProcess::BuildPhysicsTable(pd);
  BuildAsymmetryTable(pd);  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusPolarizedAnnihilation::PreparePhysicsTable(const G4ParticleDefinition& pd)
{
  G4VEmProcess::PreparePhysicsTable(pd);
  theAsymmetryTable = G4PhysicsTableHelper::PreparePhysicsTable(theAsymmetryTable);
  theTransverseAsymmetryTable = G4PhysicsTableHelper::PreparePhysicsTable(theTransverseAsymmetryTable);
}

void G4eplusPolarizedAnnihilation::BuildAsymmetryTable(const G4ParticleDefinition& part)
{
  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  G4cout<<" annih-numOfCouples="<<numOfCouples<<"\n";
  for(size_t i=0; i<numOfCouples; ++i) {
    G4cout<<"annih- "<<i<<"/"<<numOfCouples<<"\n";
    if (!theAsymmetryTable) break;
    G4cout<<"annih- "<<theAsymmetryTable->GetFlag(i)<<"\n";
    if (theAsymmetryTable->GetFlag(i)) {
     G4cout<<" building pol-annih ... \n";

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);

      // use same parameters as for lambda
      G4PhysicsVector* aVector = LambdaPhysicsVector(couple);
      G4PhysicsVector* tVector = LambdaPhysicsVector(couple);

      for (G4int j = 0 ; j < LambdaBinning() ; ++j ) {
	G4double lowEdgeEnergy = aVector->GetLowEdgeEnergy(j);
	G4double tasm=0.;
	G4double asym = ComputeAsymmetry(lowEdgeEnergy, couple, part, 0., tasm);
	aVector->PutValue(j,asym);
	tVector->PutValue(j,tasm);
      }

      G4PhysicsTableHelper::SetPhysicsVector(theAsymmetryTable, i, aVector);
      G4PhysicsTableHelper::SetPhysicsVector(theTransverseAsymmetryTable, i, tVector);
    }
  }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eplusPolarizedAnnihilation::ComputeAsymmetry(G4double energy,
			    const G4MaterialCutsCouple* couple,
			    const G4ParticleDefinition& aParticle,
			    G4double cut,
			    G4double &tAsymmetry)
{
 G4double lAsymmetry = 0.0;
 	  tAsymmetry = 0.0;

 // calculate polarized cross section
 theTargetPolarization=G4ThreeVector(0.,0.,1.);
 emModel->SetTargetPolarization(theTargetPolarization);
 emModel->SetBeamPolarization(theTargetPolarization);
 G4double sigma2=emModel->CrossSection(couple,&aParticle,energy,cut,energy);

 // calculate transversely polarized cross section
 theTargetPolarization=G4ThreeVector(1.,0.,0.);
 emModel->SetTargetPolarization(theTargetPolarization);
 emModel->SetBeamPolarization(theTargetPolarization);
 G4double sigma3=emModel->CrossSection(couple,&aParticle,energy,cut,energy);

 // calculate unpolarized cross section
 theTargetPolarization=G4ThreeVector();
 emModel->SetTargetPolarization(theTargetPolarization);
 emModel->SetBeamPolarization(theTargetPolarization);
 G4double sigma0=emModel->CrossSection(couple,&aParticle,energy,cut,energy);

 // determine assymmetries
  if (sigma0>0.) {
    lAsymmetry=sigma2/sigma0-1.;
    tAsymmetry=sigma3/sigma0-1.;
                 }
 return lAsymmetry;

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusPolarizedAnnihilation::PrintInfo()
{
  G4cout << "      Polarized model for annihilation into 2 photons"
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4eplusPolarizedAnnihilation::AtRestDoIt(const G4Track& aTrack,
                                                     const G4Step& )
//
// Performs the e+ e- annihilation when both particles are assumed at rest.
// It generates two back to back photons with energy = electron_mass.
// The angular distribution is isotropic.
// GEANT4 internal units
//
// Note : Effects due to binding of atomic electrons are negliged.
{
  fParticleChange.InitializeForPostStep(aTrack);

  fParticleChange.SetNumberOfSecondaries(2);

  G4double cosTeta = 2.*G4UniformRand()-1. , sinTeta = std::sqrt(1.-cosTeta*cosTeta);
  G4double phi     = twopi * G4UniformRand();
  G4ThreeVector direction (sinTeta*std::cos(phi), sinTeta*std::sin(phi), cosTeta);
  fParticleChange.AddSecondary( new G4DynamicParticle (G4Gamma::Gamma(),
                                            direction, electron_mass_c2) );
  fParticleChange.AddSecondary( new G4DynamicParticle (G4Gamma::Gamma(),
                                           -direction, electron_mass_c2) );
  // Kill the incident positron
  //
  fParticleChange.ProposeTrackStatus(fStopAndKill);
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
