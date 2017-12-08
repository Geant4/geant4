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
// $Id: G4ePolarizedIonisation.cc 105740 2017-08-16 13:05:44Z gcosmo $
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4ePolarizedIonisation
//
// Author:        A.Schaelicke on base of Vladimir Ivanchenko code
//
// Creation date: 10.11.2005
//
// Modifications:
//
// 10-11-05, include polarization description (A.Schaelicke)
// , create asymmetry table and determine interactionlength 
// , update polarized differential cross section 
//
// 20-08-06, modified interface (A.Schaelicke)
// 11-06-07, add PostStepGetPhysicalInteractionLength (A.Schalicke)
//
// Class Description:
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ePolarizedIonisation.hh"
#include "G4Electron.hh"
#include "G4UniversalFluctuation.hh"
#include "G4UnitsTable.hh"

#include "G4PolarizedMollerBhabhaModel.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PolarizationManager.hh"
#include "G4PolarizationHelper.hh"
#include "G4StokesVector.hh"
#include "G4EmParameters.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ePolarizedIonisation::G4ePolarizedIonisation(const G4String& name)
  : G4VEnergyLossProcess(name),
    theElectron(G4Electron::Electron()),
    isElectron(true),
    isInitialised(false),
    theTargetPolarization(0.,0.,0.),
    theAsymmetryTable(nullptr),
    theTransverseAsymmetryTable(nullptr)
{
  verboseLevel=0;
  SetProcessSubType(fIonisation);
  SetSecondaryParticle(theElectron);
  flucModel = nullptr;
  emModel = nullptr; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ePolarizedIonisation::~G4ePolarizedIonisation()
{
  CleanTables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ePolarizedIonisation::CleanTables()
{
  if(theAsymmetryTable) {
    theAsymmetryTable->clearAndDestroy();
    delete theAsymmetryTable;
    theAsymmetryTable = nullptr;
  }
  if(theTransverseAsymmetryTable) {
    theTransverseAsymmetryTable->clearAndDestroy();
    delete theTransverseAsymmetryTable;
    theTransverseAsymmetryTable = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4ePolarizedIonisation::MinPrimaryEnergy(const G4ParticleDefinition*,
					 const G4Material*, G4double cut)
{
  G4double x = cut;
  if(isElectron) { x += cut; }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4ePolarizedIonisation::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Electron::Electron() || &p == G4Positron::Positron());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ePolarizedIonisation::InitialiseEnergyLossProcess(
		    const G4ParticleDefinition* part,
		    const G4ParticleDefinition* /*part2*/)
{
  if(!isInitialised) {

    if(part == G4Positron::Positron()) { isElectron = false; }

    if (!FluctModel()) { SetFluctModel(new G4UniversalFluctuation()); }
    flucModel = FluctModel();

    emModel = new  G4PolarizedMollerBhabhaModel();
    SetEmModel(emModel);
    G4EmParameters* param = G4EmParameters::Instance();
    emModel->SetLowEnergyLimit(param->MinKinEnergy());
    emModel->SetHighEnergyLimit(param->MaxKinEnergy());
    AddEmModel(1, emModel, flucModel);

    isInitialised = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ePolarizedIonisation::PrintInfo()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ePolarizedIonisation::GetMeanFreePath(const G4Track& track,
						 G4double step,
						 G4ForceCondition* cond)
{
  // *** get unploarised mean free path from lambda table ***
  G4double mfp = G4VEnergyLossProcess::GetMeanFreePath(track, step, cond);
  if(theAsymmetryTable && theTransverseAsymmetryTable && mfp < DBL_MAX) {
    mfp *= ComputeSaturationFactor(track);
  }
  if (verboseLevel>=2) {
    G4cout << "G4ePolarizedIonisation::MeanFreePath:  " 
	   << mfp / mm << " mm " << G4endl;
  }
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ePolarizedIonisation::PostStepGetPhysicalInteractionLength(const G4Track& track,
                                              G4double step,
                                              G4ForceCondition* cond)
{
  // save previous values
  G4double nLength = theNumberOfInteractionLengthLeft;
  G4double iLength = currentInteractionLength;

  // *** get unpolarised mean free path from lambda table ***
  // this changes theNumberOfInteractionLengthLeft and currentInteractionLength
  G4double x = G4VEnergyLossProcess::PostStepGetPhysicalInteractionLength(track, step, cond);
  G4double x0 = x;
  G4double satFact = 1.;
  
  // *** add corrections on polarisation ***
  if(theAsymmetryTable && theTransverseAsymmetryTable && x < DBL_MAX) {
    satFact = ComputeSaturationFactor(track);
    G4double curLength = currentInteractionLength*satFact;
    G4double prvLength = iLength*satFact;
    if(nLength > 0.0) {
      theNumberOfInteractionLengthLeft = 
        std::max(nLength - step/prvLength, 0.0);
    }
    x = theNumberOfInteractionLengthLeft * curLength;
  }
  if (verboseLevel>=2) {
    G4cout << "G4ePolarizedIonisation::PostStepGPIL: " 
           << std::setprecision(8) << x/mm  << " mm;" << G4endl
           << "                   unpolarized value: "
           << std::setprecision(8) << x0/mm << " mm." << G4endl;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4ePolarizedIonisation::ComputeSaturationFactor(const G4Track& track)
{
  G4Material*         aMaterial = track.GetMaterial();
  G4VPhysicalVolume*  aPVolume  = track.GetVolume();
  G4LogicalVolume*    aLVolume  = aPVolume->GetLogicalVolume();
    
  G4PolarizationManager * polarizationManger = G4PolarizationManager::GetInstance();
    
  const G4bool volumeIsPolarized = polarizationManger->IsPolarized(aLVolume);
  G4StokesVector volPolarization = polarizationManger->GetVolumePolarization(aLVolume);

  G4double factor = 1.0;

  if (volumeIsPolarized && !volPolarization.IsZero()) {
     
    // *** get asymmetry, if target is polarized ***
    const G4DynamicParticle* aDynamicPart = track.GetDynamicParticle();
    const G4double energy = aDynamicPart->GetKineticEnergy();
    const G4StokesVector polarization = track.GetPolarization();
    const G4ParticleMomentum direction0 = aDynamicPart->GetMomentumDirection();

    if (verboseLevel>=2) {
      G4cout << "G4ePolarizedIonisation::ComputeSaturationFactor: " << G4endl;      
      G4cout << " Energy(MeV)  " << energy/MeV  << G4endl;
      G4cout << " Direction    " << direction0  << G4endl;
      G4cout << " Polarization " << polarization  << G4endl;
      G4cout << " MaterialPol. " << volPolarization  << G4endl;
      G4cout << " Phys. Volume " << aPVolume->GetName() << G4endl;
      G4cout << " Log. Volume  " << aLVolume->GetName() << G4endl;
      G4cout << " Material     " << aMaterial          << G4endl;
    }
    
    size_t midx = CurrentMaterialCutsCoupleIndex();
    const G4PhysicsVector* aVector = nullptr;
    const G4PhysicsVector* bVector = nullptr;
    if(midx < theAsymmetryTable->size()) { 
      aVector = (*theAsymmetryTable)(midx);
    }
    if(midx < theTransverseAsymmetryTable->size()) { 
      bVector = (*theTransverseAsymmetryTable)(midx);
    }
    if(aVector && bVector) {
      G4double lAsymmetry = aVector->Value(energy);
      G4double tAsymmetry = bVector->Value(energy);
      G4double polZZ = polarization.z()*(volPolarization*direction0);
      G4double polXX = polarization.x()*
	(volPolarization*G4PolarizationHelper::GetParticleFrameX(direction0));
      G4double polYY = polarization.y()*
	(volPolarization*G4PolarizationHelper::GetParticleFrameY(direction0));

      factor /= (1. + polZZ*lAsymmetry + (polXX + polYY)*tAsymmetry);

      if (verboseLevel>=2) {
	G4cout << " Asymmetry:     " << lAsymmetry << ", " << tAsymmetry  << G4endl;
	G4cout << " PolProduct:    " << polXX << ", " << polYY << ", " << polZZ << G4endl;
	G4cout << " Factor:        " << factor         << G4endl;
      }
    } else {
      G4ExceptionDescription ed;
      ed << "Problem with asymmetry tables: material index " << midx 
	 << " is out of range or tables are not filled";
      G4Exception("G4ePolarizedIonisation::ComputeSaturationFactor","em0048",
		  JustWarning, ed, "");
    }
  }
  return factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ePolarizedIonisation::BuildPhysicsTable(
        const G4ParticleDefinition& part)
{
  // *** build DEDX and (unpolarized) cross section tables
  G4VEnergyLossProcess::BuildPhysicsTable(part);
  G4bool master = true;
  const G4ePolarizedIonisation* masterProcess = 
    static_cast<const G4ePolarizedIonisation*>(GetMasterProcess());
  if(masterProcess && masterProcess != this) { master = false; }
  if(master) { BuildAsymmetryTables(part); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ePolarizedIonisation::BuildAsymmetryTables(
	const G4ParticleDefinition& part)
{
  // cleanup old, initialise new table
  CleanTables();
  theAsymmetryTable = 
    G4PhysicsTableHelper::PreparePhysicsTable(theAsymmetryTable);
  theTransverseAsymmetryTable = 
    G4PhysicsTableHelper::PreparePhysicsTable(theTransverseAsymmetryTable);

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  for (size_t j=0 ; j < numOfCouples; j++ ) {
    // get cut value
    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(j);

    G4double cut = (*theCoupleTable->GetEnergyCutsVector(1))[j];

    //create physics vectors then fill it (same parameters as lambda vector)
    G4PhysicsVector * ptrVectorA = LambdaPhysicsVector(couple,cut);
    G4PhysicsVector * ptrVectorB = LambdaPhysicsVector(couple,cut);
    size_t bins = ptrVectorA->GetVectorLength();

    for (size_t i = 0 ; i < bins ; i++ ) {
      G4double lowEdgeEnergy = ptrVectorA->Energy(i);
      G4double tasm=0.;
      G4double asym = ComputeAsymmetry(lowEdgeEnergy, couple, part, cut, tasm);
      ptrVectorA->PutValue(i,asym);
      ptrVectorB->PutValue(i,tasm);
    }
    theAsymmetryTable->insertAt( j , ptrVectorA ) ;
    theTransverseAsymmetryTable->insertAt( j , ptrVectorB ) ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4ePolarizedIonisation::ComputeAsymmetry(G4double energy,
					 const G4MaterialCutsCouple* couple,
					 const G4ParticleDefinition& aParticle,
					 G4double cut,
					 G4double & tAsymmetry)
{
  G4double lAsymmetry = 0.0;
  	   tAsymmetry = 0.0;
  if (isElectron) { lAsymmetry = tAsymmetry = -1.0; }

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
  G4double sigma0 = emModel->CrossSection(couple,&aParticle,energy,cut,energy);
  // determine assymmetries
  if (sigma0 > 0.) {
    lAsymmetry=sigma2/sigma0 - 1.;
    tAsymmetry=sigma3/sigma0 - 1.;
  }
  if (std::fabs(lAsymmetry)>1.) {
    G4cout<<"G4ePolarizedIonisation::ComputeAsymmetry WARNING: E(MeV)= " 
	  << energy << " lAsymmetry= "<<lAsymmetry
	  <<" ("<<std::fabs(lAsymmetry)-1.<<")\n";
  }
  if (std::fabs(tAsymmetry)>1.) {
    G4cout<<" energy="<<energy<<"\n";
    G4cout<<"G4ePolarizedIonisation::ComputeAsymmetry WARNING: E(MeV)= " 
	  << energy << " tAsymmetry= "<<tAsymmetry
	  <<" ("<<std::fabs(tAsymmetry)-1.<<")\n";
  }
  return lAsymmetry;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

