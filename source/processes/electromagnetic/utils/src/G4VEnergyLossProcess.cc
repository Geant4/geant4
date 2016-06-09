//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4VEnergyLossProcess.cc,v 1.57 2005/05/27 18:38:33 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4VEnergyLossProcess
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.01.2002
//
// Modifications:
//
// 13-11-02 Minor fix - use normalised direction (V.Ivanchenko)
// 04-12-02 Minor change in PostStepDoIt (V.Ivanchenko)
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 26-12-02 Secondary production moved to derived classes (V.Ivanchenko)
// 04-01-03 Fix problem of very small steps for ions (V.Ivanchenko)
// 20-01-03 Migrade to cut per region (V.Ivanchenko)
// 24-01-03 Temporarily close a control on usage of couples (V.Ivanchenko)
// 24-01-03 Make models region aware (V.Ivanchenko)
// 05-02-03 Fix compilation warnings (V.Ivanchenko)
// 06-02-03 Add control on tmax in PostStepDoIt (V.Ivanchenko)
// 13-02-03 SubCutoffProcessors defined for regions (V.Ivanchenko)
// 15-02-03 Lambda table can be scaled (V.Ivanchenko)
// 17-02-03 Fix problem of store/restore tables (V.Ivanchenko)
// 18-02-03 Add control on CutCouple usage (V.Ivanchenko)
// 26-02-03 Simplify control on GenericIons (V.Ivanchenko)
// 06-03-03 Control on GenericIons using SubType + update verbose (V.Ivanchenko)
// 10-03-03 Add Ion registration (V.Ivanchenko)
// 22-03-03 Add Initialisation of cash (V.Ivanchenko)
// 26-03-03 Remove finalRange modification (V.Ivanchenko)
// 09-04-03 Fix problem of negative range limit for non integral (V.Ivanchenko)
// 26-04-03 Fix retrieve tables (V.Ivanchenko)
// 06-05-03 Set defalt finalRange = 1 mm (V.Ivanchenko)
// 12-05-03 Update range calculations + lowKinEnergy (V.Ivanchenko)
// 13-05-03 Add calculation of precise range (V.Ivanchenko)
// 23-05-03 Remove tracking cuts (V.Ivanchenko)
// 03-06-03 Fix initialisation problem for STD ionisation (V.Ivanchenko)
// 21-07-03 Add UpdateEmModel method (V.Ivanchenko)
// 03-11-03 Fix initialisation problem in RetrievePhysicsTable (V.Ivanchenko)
// 04-11-03 Add checks in RetrievePhysicsTable (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 21-01-04 Migrade to G4ParticleChangeForLoss (V.Ivanchenko)
// 27-02-04 Fix problem of loss in low presure gases, cleanup precise range
//          calculation, use functions ForLoss in AlongStepDoIt (V.Ivanchenko)
// 10-03-04 Fix a problem of Precise Range table (V.Ivanchenko)
// 19-03-04 Fix a problem energy below lowestKinEnergy (V.Ivanchenko)
// 31-03-04 Fix a problem of retrieve tables (V.Ivanchenko)
// 21-07-04 Check weather AtRest are active or not (V.Ivanchenko)
// 03-08-04 Add pointer of DEDX table to all processes (V.Ivanchenko)
// 06-08-04 Clear up names of member functions (V.Ivanchenko)
// 06-08-04 Clear up names of member functions (V.Ivanchenko)
// 27-08-04 Add NeedBuildTables method (V.Ivanchneko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 11-03-05 Shift verbose level by 1 (V.Ivantchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 11-04-05 Use MaxSecondaryEnergy from a model (V.Ivanchenko)
//
// Class Description:
//
// It is the unified energy loss process it calculates the continuous
// energy loss for charged particles using a set of Energy Loss
// models valid for different energy regions. There are a possibility
// to create and access to dE/dx and range tables, or to calculate
// that information on fly.
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VEnergyLossProcess.hh"
#include "G4LossTableManager.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4DataVector.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4VParticleChange.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4VSubCutoffProcessor.hh"
#include "G4ProcessManager.hh"
#include "G4UnitsTable.hh"
#include "G4GenericIon.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysicsTableHelper.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossProcess::G4VEnergyLossProcess(const G4String& name, G4ProcessType type):
                 G4VContinuousDiscreteProcess(name, type),
  nSCoffRegions(0),
  idxSCoffRegions(0),
  theDEDXTable(0),
  theRangeTableForLoss(0),
  theDEDXunRestrictedTable(0),
  thePreciseRangeTable(0),
  theSecondaryRangeTable(0),
  theInverseRangeTable(0),
  theLambdaTable(0),
  theSubLambdaTable(0),
  theDEDXAtMaxEnergy(0),
  theRangeAtMaxEnergy(0),
  theEnergyOfCrossSectionMax(0),
  theCrossSectionMax(0),
  particle(0),
  baseParticle(0),
  secondaryParticle(0),
  currentCouple(0),
  nDEDXBins(90),
  nDEDXBinsForRange(70),
  nLambdaBins(90),
  linLossLimit(0.05),
  minSubRange(0.1),
  defaultRoverRange(0.2),
  defaultIntegralRange(1.0),
  lambdaFactor(0.1),
  mfpKinEnergy(0.0),
  lossFluctuationFlag(true),
  lossFluctuationArePossible(true),
  rndmStepFlag(false),
  tablesAreBuilt(false),
  integral(true),
  meanFreePath(false),
  aboveCSmax(true),
  isIonisation(true)
{

  lowestKinEnergy      = 1.*eV;
  minKinEnergy         = 0.1*keV;
  maxKinEnergy         = 100.0*GeV;
  maxKinEnergyForRange = 1.0*GeV;

  pParticleChange = &fParticleChange;

  // default dRoverRange and finalRange
  SetStepFunction(defaultIntegralRange, 1.0*mm);
  SetVerboseLevel(1);

  modelManager = new G4EmModelManager();
  (G4LossTableManager::Instance())->Register(this);
  scoffProcessors.clear();
  scoffRegions.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossProcess::~G4VEnergyLossProcess()
{
  Clear();

  if ( !baseParticle ) {
    if(theDEDXTable && theRangeTableForLoss) theDEDXTable->clearAndDestroy();
    if(theDEDXunRestrictedTable && thePreciseRangeTable)
       theDEDXunRestrictedTable->clearAndDestroy();
    if(thePreciseRangeTable) thePreciseRangeTable->clearAndDestroy();
    if(theRangeTableForLoss) theRangeTableForLoss->clearAndDestroy();
    if(theInverseRangeTable) theInverseRangeTable->clearAndDestroy();
    if(theLambdaTable) theLambdaTable->clearAndDestroy();
    if(theSubLambdaTable) theSubLambdaTable->clearAndDestroy();
  }

  if (nSCoffRegions) {
    for (G4int i=0; i<nSCoffRegions; i++) {
      if (scoffProcessors[i]) {
	for (G4int j=i+1; j<nSCoffRegions; j++) {
	  if(scoffProcessors[i] == scoffProcessors[j]) scoffProcessors[j] = 0;
	}
        delete scoffProcessors[i];
      }
    }
    scoffProcessors.clear();
    scoffRegions.clear();
  }
  delete modelManager;
  (G4LossTableManager::Instance())->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::Clear()
{
  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::Clear() for " << GetProcessName() << G4endl;
  }

  if(theDEDXAtMaxEnergy) delete [] theDEDXAtMaxEnergy;
  if(theRangeAtMaxEnergy) delete [] theRangeAtMaxEnergy;
  if(theEnergyOfCrossSectionMax) delete [] theEnergyOfCrossSectionMax;
  if(theCrossSectionMax) delete [] theCrossSectionMax;

  theDEDXAtMaxEnergy = 0;
  theRangeAtMaxEnergy = 0;
  theEnergyOfCrossSectionMax = 0,
  theCrossSectionMax = 0,
  tablesAreBuilt = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::PreparePhysicsTable(const G4ParticleDefinition& part)
{

  // Are particle defined?
  if( !particle ) {
    if(part.GetParticleType() == "nucleus" && part.GetParticleSubType() == "generic") 
         particle = G4GenericIon::GenericIon();
    else particle = &part;
  }

  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::PreparePhysicsTable for "
           << GetProcessName()
           << " for " << part.GetParticleName()
           << " local: " << particle->GetParticleName()
           << G4endl;
  }

  G4LossTableManager* lManager = G4LossTableManager::Instance();

  if (&part != particle) {
    if (part.GetParticleType() == "nucleus") lManager->RegisterIon(&part, this);
    else                                     lManager->RegisterExtraParticle(&part, this);
    return;
  }

  Clear();

  currentCouple = 0;
  preStepLambda = 0.0;
  mfpKinEnergy  = DBL_MAX;
  preStepMFP    = DBL_MAX;

  // Base particle and set of models can be defined here
  InitialiseEnergyLossProcess(particle, baseParticle);

  // Tables preparation
  if (!baseParticle) {
    
    theDEDXTable = G4PhysicsTableHelper::PreparePhysicsTable(theDEDXTable);
    if (lManager->BuildPreciseRange()) {
      theDEDXunRestrictedTable = G4PhysicsTableHelper::PreparePhysicsTable(theDEDXunRestrictedTable);
      if (isIonisation)
        thePreciseRangeTable = G4PhysicsTableHelper::PreparePhysicsTable(thePreciseRangeTable);
    }

    if (isIonisation) {
      theRangeTableForLoss = G4PhysicsTableHelper::PreparePhysicsTable(theRangeTableForLoss);
      theInverseRangeTable = G4PhysicsTableHelper::PreparePhysicsTable(theInverseRangeTable);
    }
    theLambdaTable       = G4PhysicsTableHelper::PreparePhysicsTable(theLambdaTable);
    if (nSCoffRegions)
      theSubLambdaTable = G4PhysicsTableHelper::PreparePhysicsTable(theSubLambdaTable);
  }

  G4double initialCharge = particle->GetPDGCharge();
  G4double initialMass   = particle->GetPDGMass();
  chargeSquare = initialCharge*initialCharge/(eplus*eplus);
  chargeSqRatio = 1.0;
  massRatio = 1.0;
  reduceFactor = 1.0;

  if (baseParticle) {
    massRatio = (baseParticle->GetPDGMass())/initialMass;
    G4double q = initialCharge/baseParticle->GetPDGCharge();
    chargeSqRatio = q*q;
    reduceFactor = 1.0/(chargeSqRatio*massRatio);
  }

  theCuts = modelManager->Initialise(particle, secondaryParticle, minSubRange, verboseLevel);

  // Sub Cutoff Regime

  idxSCoffRegions.clear();

  const G4ProductionCutsTable* theCoupleTable=
          G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  if (nSCoffRegions) {
    const G4DataVector* theSubCuts = modelManager->SubCutoff();

    for (G4int i=0; i<nSCoffRegions; i++) {
      scoffProcessors[i]->Initialise(particle, secondaryParticle, theCuts, theSubCuts);
    }
    for (size_t j=0; j<numOfCouples; j++) {

      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(j);
      const G4ProductionCuts* pcuts = couple->GetProductionCuts();
      G4int reg = nSCoffRegions;
      do {reg--;} while (reg && pcuts != (scoffRegions[reg]->GetProductionCuts()));
      idxSCoffRegions.push_back(reg);
    }
  }

  lManager->EnergyLossProcessIsInitialised(particle, this);

  if (1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::Initialise() is done "
           << " chargeSqRatio= " << chargeSqRatio
           << " massRatio= " << massRatio
           << " reduceFactor= " << reduceFactor << G4endl;
    if (nSCoffRegions) {
      G4cout << " SubCutoff Regime is ON for regions: " << G4endl;
      for (G4int i=0; i<nSCoffRegions; i++) {
        const G4Region* r = scoffRegions[i];
	G4cout << "           " << r->GetName() << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  if(1 < verboseLevel) {
  //    G4cout << "========================================================" << G4endl;
    G4cout << "### G4VEnergyLossProcess::BuildPhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << "; local: " << particle->GetParticleName();
    if(baseParticle) G4cout << "; base: " << baseParticle->GetParticleName();
    G4cout << G4endl;
  }

  if(!tablesAreBuilt && &part == particle)
    G4LossTableManager::Instance()->BuildPhysicsTable(particle, this);

  if(0 < verboseLevel && (&part == particle) && !baseParticle) PrintInfoDefinition();

  if(1 < verboseLevel) {
    G4cout << "### G4VEnergyLossProcess::BuildPhysicsTable() done for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::AddEmModel(G4int order, G4VEmModel* p, G4VEmFluctuationModel* fluc,
                                const G4Region* region)
{
  modelManager->AddEmModel(order, p, fluc, region);
  if(p) p->SetParticleChange(pParticleChange, fluc);
  if(!fluc) {
    lossFluctuationFlag = false;
    lossFluctuationArePossible = false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::UpdateEmModel(const G4String& nam, G4double emin, G4double emax)
{
  modelManager->UpdateEmModel(nam, emin, emax);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::AddSubCutoffProcessor(G4VSubCutoffProcessor* p,
                                           const G4Region* r)
{
  if( !p ) {
    G4cout << "G4VEnergyLossProcess::AddSubCutoffProcessor WARNING: no SubCutoffProcessor defined."
           << G4endl;
    return;
  }
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  if (!r) r = regionStore->GetRegion("DefaultRegionForTheWorld", false);
  if (nSCoffRegions) {
    for (G4int i=0; i<nSCoffRegions; i++) {
      if (r == scoffRegions[i]) {
        if ( scoffProcessors[i] ) delete scoffProcessors[i];
	scoffProcessors[i] = p;
        return;
      }
    }
  }
  scoffProcessors.push_back(p);
  scoffRegions.push_back(r);
  nSCoffRegions++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossProcess::BuildDEDXTable()
{

  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTable() for "
           << GetProcessName()
           << " and particle " << particle->GetParticleName()
           << G4endl;
  }

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  if(1 < verboseLevel) {
    G4cout << numOfCouples << " materials"
           << " minKinEnergy= " << minKinEnergy
           << " maxKinEnergy= " << maxKinEnergy
           << G4endl;
  }

  for(size_t i=0; i<numOfCouples; i++) {

    if(2 < verboseLevel) 
      G4cout << "G4VEnergyLossProcess::BuildDEDXVector flag=  " << theDEDXTable->GetFlag(i) << G4endl;

    if (theDEDXTable->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
      G4PhysicsVector* aVector = DEDXPhysicsVector(couple);
      modelManager->FillDEDXVector(aVector, couple);

      // Insert vector for this material into the table
      G4PhysicsTableHelper::SetPhysicsVector(theDEDXTable, i, aVector);
    }
  }

  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTable(): table is built for "
           << particle->GetParticleName()
           << G4endl;
    if(2 < verboseLevel) G4cout << (*theDEDXTable) << G4endl;
  }

  return theDEDXTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossProcess::BuildDEDXTableForPreciseRange()
{

  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTableForPreciseRange() for "
           << GetProcessName()
           << " and particle " << particle->GetParticleName()
           << G4endl;
  }

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  if(0 < verboseLevel) {
    G4cout << numOfCouples << " materials"
           << " minKinEnergy= " << minKinEnergy
           << " maxKinEnergy= " << maxKinEnergy
           << G4endl;
  }

  for(size_t i=0; i<numOfCouples; i++) {

    if (theDEDXunRestrictedTable->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
      G4PhysicsVector* aVector = DEDXPhysicsVectorForPreciseRange(couple);
      modelManager->FillDEDXVectorForPreciseRange(aVector, couple);

      // Insert vector for this material into the table
      G4PhysicsTableHelper::SetPhysicsVector(theDEDXunRestrictedTable, i, aVector);
    }
  }

  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTableForPreciseRange(): table is built for "
           << particle->GetParticleName()
           << G4endl;
    if(2 < verboseLevel) G4cout << (*theDEDXunRestrictedTable) << G4endl;
  }

  return theDEDXunRestrictedTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossProcess::BuildLambdaTable()
{

  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildLambdaTable() for process "
           << GetProcessName() << " and particle "
           << particle->GetParticleName()
           << G4endl;
  }

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  for(size_t i=0; i<numOfCouples; i++) {

    if (theLambdaTable->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
      G4PhysicsVector* aVector = LambdaPhysicsVector(couple);
      modelManager->FillLambdaVector(aVector, couple);

      // Insert vector for this material into the table
      G4PhysicsTableHelper::SetPhysicsVector(theLambdaTable, i, aVector);
    }
  }

  if(1 < verboseLevel) {
    G4cout << "Lambda table is built for "
           << particle->GetParticleName()
           << G4endl;
  }

  return theLambdaTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossProcess::BuildLambdaSubTable()
{
  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildLambdaSubTable() for process "
           << GetProcessName() << " and particle "
           << particle->GetParticleName() << G4endl;
  }

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  for(size_t i=0; i<numOfCouples; i++) {

    if (theSubLambdaTable->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
      G4PhysicsVector* aVector = SubLambdaPhysicsVector(couple);
      modelManager->FillSubLambdaVector(aVector, couple);

      // Insert vector for this material into the table
      G4PhysicsTableHelper::SetPhysicsVector(theSubLambdaTable, i, aVector);
    }
  }

  if(1 < verboseLevel) {
    G4cout << "Table is built for "
           << particle->GetParticleName()
           << G4endl;
  }

  return theSubLambdaTable;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEnergyLossProcess::AlongStepDoIt(const G4Track& track,
                                                       const G4Step& step)
{
  fParticleChange.InitializeForAlongStep(track);
  // The process has range table - calculate energy loss
  if(!theRangeTableForLoss) return &fParticleChange;

  // Get the actual (true) Step length
  G4double length = step.GetStepLength();
  G4double eloss  = 0.0;

  /*
  if(-1 < verboseLevel) {
    const G4ParticleDefinition* d = track.GetDefinition();
    G4cout << "AlongStepDoIt for "
           << GetProcessName() << " and particle "
           << d->GetParticleName()
           << "  eScaled(MeV)= " << preStepScaledEnergy/MeV
           << "  slim(mm)= " << fRange/mm
           << "  s(mm)= " << length/mm
           << "  q^2= " << chargeSqRatio
           << " md= " << d->GetPDGMass()
           << G4endl;
  }
  */

  // stopping
  if (length >= fRange) {
    eloss = preStepKinEnergy;

  // Short step
  } else if( length <= linLossLimit * fRange ) {
    eloss = GetDEDXForScaledEnergy(preStepScaledEnergy)*length;

  // Long step
  } else {
    G4double r = GetScaledRangeForScaledEnergy(preStepScaledEnergy);
    G4double x = r - length/reduceFactor;
    if(x < 0.0) {
      G4cout << "WARNING! G4VEnergyLossProcess::AlongStepDoIt: x= " << x
             << " for eScaled(MeV)= " << preStepScaledEnergy/MeV
             << " step(mm)= " << length/mm
             << " for " << track.GetDefinition()->GetParticleName()
             << G4endl;
      x = 0.0;
    }
    eloss = (ScaledKinEnergyForLoss(r) - ScaledKinEnergyForLoss(x))/massRatio;

    /*
    if(-1 < verboseLevel) {
      G4cout << "Long STEP: rPre(mm)= " << r/mm
             << " rPost(mm)= " << x/mm
             << " ePre(MeV)= " << preStepScaledEnergy/MeV
             << " eloss(MeV)= " << eloss/MeV
             << " eloss0(MeV)= " 
             << GetDEDXForScaledEnergy(preStepScaledEnergy)*length/MeV
             << G4endl;
    }
    */

  }

  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  G4VEmModel* currentModel = SelectModel(preStepScaledEnergy);
  G4double tmax = currentModel->MaxSecondaryKinEnergy(dynParticle);
  tmax = std::min(tmax,(*theCuts)[currentMaterialIndex]);

  /*  
  G4double eloss0 = eloss;
  if(-1 < verboseLevel) {
    G4cout << "Before fluct: eloss(MeV)= " << eloss/MeV
           << " tmax= " << tmax
           << " e-eloss= " << preStepKinEnergy-eloss
           << "  fluct= " << lossFluctuationFlag 
           << G4endl;
  }
  */

  // Sample fluctuations
  if (lossFluctuationFlag && eloss < preStepKinEnergy) {

    eloss = currentModel->GetModelOfFluctuations()->
      SampleFluctuations(currentMaterial,dynParticle,tmax,length,eloss);
  }
  /*
  if(-1 < verboseLevel) {
    G4cout << "After fluct: eloss(MeV)= " << eloss/MeV
           << " fluc= " << (eloss-eloss0)/MeV
           << " currentChargeSquare= " << chargeSquare
           << " massRatio= " << massRatio
           << G4endl;
  }
  */

  // Corrections, which cannot be tabulated 
  CorrectionsAlongStep(currentCouple, dynParticle, eloss, length);
 
  G4double finalT = preStepKinEnergy - eloss;
  if (finalT <= lowestKinEnergy) finalT = 0.0;
  eloss = preStepKinEnergy-finalT;

  fParticleChange.SetProposedKineticEnergy(finalT);
  fParticleChange.ProposeLocalEnergyDeposit(eloss);

  /*  
  if(-1 < verboseLevel) {
    G4cout << "Final value eloss(MeV)= " << eloss/MeV
           << " preStepKinEnergy= " << preStepKinEnergy
           << " postStepKinEnergy= " << finalT
           << " lossFlag= " << lossFluctuationFlag
           << G4endl;
  }
  */

  return &fParticleChange;
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEnergyLossProcess::PostStepDoIt(const G4Track& track,
                                                  const G4Step& step)
{
  fParticleChange.InitializeForPostStep(track);
  G4double finalT = track.GetKineticEnergy();
  G4double postStepScaledEnergy = finalT*massRatio;

  // Integral approach
  if (integral) {
    G4double lx = GetLambdaForScaledEnergy(postStepScaledEnergy);
    if(preStepLambda<lx && 1 < verboseLevel) {
      G4cout << "WARING: for " << particle->GetParticleName()
             << " and " << GetProcessName()
             << " E(MeV)= " << finalT/MeV
             << " preLambda= " << preStepLambda << " < " << lx << " (postLambda) "
	     << G4endl;
    }
    if(preStepLambda*G4UniformRand() > lx)
      return G4VContinuousDiscreteProcess::PostStepDoIt(track,step);
  }

  G4VEmModel* currentModel = SelectModel(postStepScaledEnergy);
  G4double tmax = (*theCuts)[currentMaterialIndex];

  std::vector<G4DynamicParticle*>* newp = SecondariesPostStep(
       currentModel, currentCouple, track.GetDynamicParticle(), tmax);

  if (newp) {
    G4int num = newp->size();
    fParticleChange.SetNumberOfSecondaries(num);
    for (G4int i=0; i<num; i++) {
      fParticleChange.AddSecondary((*newp)[i]);
    }
    delete newp;
  }

  /*
  if(-1 < verboseLevel) {
    const G4ParticleDefinition* pd = track.GetDynamicParticle()->GetDefinition();
    G4cout << GetProcessName()
           << "::PostStepDoIt: Sample secondary; E= " << finalT/MeV
           << " MeV; model= (" << currentModel->LowEnergyLimit(pd)
           << ", " <<  currentModel->HighEnergyLimit(pd) << ")"
           << "  preStepLambda= " << preStepLambda
           << G4endl;
  }
  */
  
  finalT = fParticleChange.GetProposedKineticEnergy();
  if (finalT <= lowestKinEnergy) {
    fParticleChange.SetProposedKineticEnergy(0.0);
    return &fParticleChange;
  }

  return G4VContinuousDiscreteProcess::PostStepDoIt(track,step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::PrintInfoDefinition()
{
  if(0 < verboseLevel) {
    G4cout << G4endl << GetProcessName() << ":   tables are built for  "
           << particle->GetParticleName()
           << G4endl
           << "      dE/dx and range tables from "
 	   << G4BestUnit(minKinEnergy,"Energy")
           << " to " << G4BestUnit(maxKinEnergy,"Energy")
           << " in " << nDEDXBins << " bins." << G4endl
           << "      Lambda tables from threshold to "
           << G4BestUnit(maxKinEnergy,"Energy")
           << " in " << nLambdaBins << " bins."
           << G4endl;
    PrintInfo();
    if(theRangeTableForLoss) {
      G4cout << "      Step function: finalRange(mm)= " << finalRange/mm
             << ", dRoverRange= " << dRoverRange
             << ", integral: " << integral
             << G4endl;
    }
    if(thePreciseRangeTable) {
      G4cout << "      Precise range table up"
             << " to " << G4BestUnit(maxKinEnergyForRange,"Energy")
             << " in " << nDEDXBinsForRange << " bins." << G4endl;
    }

    if(2 < verboseLevel) {
      G4cout << "DEDXTable address= " << theDEDXTable << G4endl;
      if(theDEDXTable) G4cout << (*theDEDXTable) << G4endl;
      G4cout << "non restricted DEDXTable address= " << theDEDXunRestrictedTable << G4endl;
      if(theDEDXunRestrictedTable) G4cout << (*theDEDXunRestrictedTable) << G4endl;
      G4cout << "PreciseRangeTable address= " << thePreciseRangeTable << G4endl;
      if(thePreciseRangeTable) G4cout << (*thePreciseRangeTable) << G4endl;
      G4cout << "RangeTableForLoss address= " << theRangeTableForLoss << G4endl;
      if(theRangeTableForLoss) G4cout << (*theRangeTableForLoss) << G4endl;
      G4cout << "InverseRangeTable address= " << theInverseRangeTable << G4endl;
      if(theInverseRangeTable) G4cout << (*theInverseRangeTable) << G4endl;
      G4cout << "LambdaTable address= " << theLambdaTable << G4endl;
      if(theLambdaTable) G4cout << (*theLambdaTable) << G4endl;
      G4cout << "SubLambdaTable address= " << theSubLambdaTable << G4endl;
      if(theSubLambdaTable) G4cout << (*theSubLambdaTable) << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetDEDXTable(G4PhysicsTable* p)
{
  if(theDEDXTable != p) theDEDXTable = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetDEDXunRestrictedTable(G4PhysicsTable* p)
{
  if(theDEDXunRestrictedTable != p) theDEDXunRestrictedTable = p;
  if(p) {
    size_t n = p->length();
    G4PhysicsVector* pv = (*p)[0];
    G4double emax = maxKinEnergyForRange;
    G4bool b;
    theDEDXAtMaxEnergy = new G4double [n];

    for (size_t i=0; i<n; i++) {
      pv = (*p)[i];
      G4double dedx = pv->GetValue(emax, b);
      theDEDXAtMaxEnergy[i] = dedx;
      //G4cout << "i= " << i << " emax(MeV)= " << emax/MeV<< " dedx= " << dedx << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetPreciseRangeTable(G4PhysicsTable* p)
{
  if(thePreciseRangeTable != p) thePreciseRangeTable = p;

  if(p) {
    size_t n = p->length();
    G4PhysicsVector* pv = (*p)[0];
    G4double emax = maxKinEnergyForRange;
    G4bool b;
    theRangeAtMaxEnergy = new G4double [n];

    for (size_t i=0; i<n; i++) {
      pv = (*p)[i];
      G4double r2 = pv->GetValue(emax, b);
      theRangeAtMaxEnergy[i] = r2;
      //G4cout << "i= " << i << " e2(MeV)= " << emax/MeV << " r2= " << r2<< G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetRangeTableForLoss(G4PhysicsTable* p)
{
  if(theRangeTableForLoss != p) {
    theRangeTableForLoss = p;
    if(1 < verboseLevel) {
      G4cout << "### Set Range table " << p << " for " << particle->GetParticleName()
             << " and process " << GetProcessName() << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetSecondaryRangeTable(G4PhysicsTable* p)
{
  theSecondaryRangeTable = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetInverseRangeTable(G4PhysicsTable* p)
{
  if(theInverseRangeTable != p) {
    theInverseRangeTable = p;
    if(1 < verboseLevel) {
      G4cout << "### Set InverseRange table " << p << " for " << particle->GetParticleName()
             << " and process " << GetProcessName() << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetLambdaTable(G4PhysicsTable* p)
{
  if(1 < verboseLevel) {
    G4cout << "### Set Lambda table " << p << " for " << particle->GetParticleName()
           << " and process " << GetProcessName() << G4endl;
  }
  if(theLambdaTable != p) theLambdaTable = p;
  tablesAreBuilt = true;

  if(p) {
    size_t n = p->length();
    G4PhysicsVector* pv = (*p)[0];
    G4double e, s, smax, emax;
    theEnergyOfCrossSectionMax = new G4double [n];
    theCrossSectionMax = new G4double [n];
    G4bool b;

    for (size_t i=0; i<n; i++) {
      pv = (*p)[i];
      emax = DBL_MAX;
      smax = 0.0;
      if(pv) {
        size_t nb = pv->GetVectorLength();
        emax = pv->GetLowEdgeEnergy(nb);
	for (size_t j=0; j<nb; j++) {
	  e = pv->GetLowEdgeEnergy(j);
	  s = pv->GetValue(e,b);
	  if(s > smax) {
	    smax = s;
	    emax = e;
	  }
	}
      }
      theEnergyOfCrossSectionMax[i] = emax;
      theCrossSectionMax[i] = smax;
      if(1 < verboseLevel) {
        G4cout << "For " << particle->GetParticleName()
               << " Max CS at i= " << i << " emax(MeV)= " << emax/MeV
               << " lambda= " << smax << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetSubLambdaTable(G4PhysicsTable* p)
{
  if(theSubLambdaTable != p) theSubLambdaTable = p;

  if (nSCoffRegions) {
    for (G4int i=0; i<nSCoffRegions; i++) {
      scoffProcessors[i]->SetLambdaSubTable(theSubLambdaTable);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEnergyLossProcess::DEDXPhysicsVector(const G4MaterialCutsCouple*)
{
  G4int nbins = nDEDXBins;
  G4PhysicsVector* v = new G4PhysicsLogVector(minKinEnergy, maxKinEnergy, nbins);
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEnergyLossProcess::DEDXPhysicsVectorForPreciseRange(
                             const G4MaterialCutsCouple*)
{
  G4int nbins = nDEDXBinsForRange;
  G4PhysicsVector* v = new G4PhysicsLogVector(minKinEnergy, maxKinEnergyForRange, nbins);
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEnergyLossProcess::LambdaPhysicsVector(const G4MaterialCutsCouple* couple)
{
  G4double cut  = (*theCuts)[couple->GetIndex()];
  G4int nbins = nLambdaBins;
  G4double tmin = std::max(MinPrimaryEnergy(particle, couple->GetMaterial(), cut),
                               minKinEnergy);
  if(tmin >= maxKinEnergy) tmin = 0.5*maxKinEnergy;
  G4PhysicsVector* v = new G4PhysicsLogVector(tmin, maxKinEnergy, nbins);
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEnergyLossProcess::SubLambdaPhysicsVector(const G4MaterialCutsCouple* couple)
{
  return LambdaPhysicsVector(couple);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::MicroscopicCrossSection(G4double kineticEnergy,
                                             const G4MaterialCutsCouple* couple)
{
  // Cross section per atom is calculated
  DefineMaterial(couple);
  G4double cross = 0.0;
  G4bool b;
  if(theLambdaTable) {
    cross = (((*theLambdaTable)[currentMaterialIndex])->
                           GetValue(kineticEnergy, b));

    cross /= currentMaterial->GetTotNbOfAtomsPerVolume();
  }

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::MeanFreePath(const G4Track& track,
                                              G4double s,
                                              G4ForceCondition* cond)
{
  return GetMeanFreePath(track, s, cond);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::ContinuousStepLimit(const G4Track& track,
                                               G4double x, G4double y, G4double& z)
{
  return GetContinuousStepLimit(track, x, y, z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetStepLimits(G4double v1, G4double v2)
{
  dRoverRange = v1;
  finalRange = v2;
  if (dRoverRange > 1.0) dRoverRange = 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetIntegral(G4bool val)
{
  if(integral != val) {
    if(val) dRoverRange = defaultIntegralRange;
    else    dRoverRange = defaultRoverRange;
  }
  integral = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetStepFunction(G4double v1, G4double v2)
{
  dRoverRange = v1;
  finalRange = v2;
  if (dRoverRange > 0.999) dRoverRange = 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetBaseParticle(const G4ParticleDefinition* p)
{
  baseParticle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetSecondaryParticle(const G4ParticleDefinition* p)
{
  secondaryParticle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4VEnergyLossProcess::StorePhysicsTable(const G4ParticleDefinition* part,
			 	               const G4String& directory,
				                     G4bool ascii)
{
  G4bool res = true;
  if ( baseParticle || part != particle ) return res;

  if ( theDEDXTable && theRangeTableForLoss ) {
    const G4String name = GetPhysicsTableFileName(part,directory,"DEDX",ascii);
    if( !theDEDXTable->StorePhysicsTable(name,ascii)) res = false;
  }

  if ( theDEDXunRestrictedTable && thePreciseRangeTable ) {
    const G4String name = GetPhysicsTableFileName(part,directory,"DEDXnr",ascii);
    if( !theDEDXTable->StorePhysicsTable(name,ascii)) res = false;
  }

  if ( thePreciseRangeTable ) {
    const G4String name = GetPhysicsTableFileName(part,directory,"PreciseRange",ascii);
    if( !thePreciseRangeTable->StorePhysicsTable(name,ascii)) res = false;
  }

  if ( theRangeTableForLoss ) {
    const G4String name = GetPhysicsTableFileName(part,directory,"Range",ascii);
    if( !theRangeTableForLoss->StorePhysicsTable(name,ascii)) res = false;
  }

  if ( theInverseRangeTable ) {
    const G4String name = GetPhysicsTableFileName(part,directory,"InverseRange",ascii);
    if( !theInverseRangeTable->StorePhysicsTable(name,ascii)) res = false;
  }

  if ( theLambdaTable ) {
    const G4String name = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
    if( !theLambdaTable->StorePhysicsTable(name,ascii)) res = false;
  }

  if ( theSubLambdaTable ) {
    const G4String name = GetPhysicsTableFileName(part,directory,"SubLambda",ascii);
    if( !theSubLambdaTable->StorePhysicsTable(name,ascii)) res = false;
  }
  if ( res ) {
    if(0 < verboseLevel) {
      G4cout << "Physics tables are stored for " << particle->GetParticleName()
             << " and process " << GetProcessName()
	     << " in the directory <" << directory
	     << "> " << G4endl;
    }
  } else {
    G4cout << "Fail to store Physics Tables for " << particle->GetParticleName()
           << " and process " << GetProcessName()
	   << " in the directory <" << directory
	   << "> " << G4endl;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4VEnergyLossProcess::RetrievePhysicsTable(const G4ParticleDefinition* part,
			  	                  const G4String& directory,
			  	                        G4bool ascii)
{
  G4bool res = true;
  const G4String particleName = part->GetParticleName();

  if(1 < verboseLevel) {
 //   G4cout << "========================================================" << G4endl;
    G4cout << "G4VEnergyLossProcess::RetrievePhysicsTable() for "
           << particleName << " and process " << GetProcessName()
           << "; tables_are_built= " << tablesAreBuilt
           << G4endl;
  }
  if(particle == part) {

    G4bool yes = true;
    G4bool fpi = true;
    if ( !baseParticle ) {
      G4String filename;

      filename = GetPhysicsTableFileName(part,directory,"DEDX",ascii);
      yes = theDEDXTable->ExistPhysicsTable(filename);
      if(yes) yes = G4PhysicsTableHelper::RetrievePhysicsTable(theDEDXTable,filename,ascii);
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "DEDX table for " << particleName << " is Retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
        fpi = false;
        if (1 < verboseLevel) {
          G4cout << "DEDX table for " << particleName << " from file <"
                 << filename << "> is not Retrieved"
                 << G4endl;
        }
      }

      filename = GetPhysicsTableFileName(part,directory,"Range",ascii);
      yes = theRangeTableForLoss->ExistPhysicsTable(filename);
      if(yes) yes = G4PhysicsTableHelper::RetrievePhysicsTable(theRangeTableForLoss,filename,ascii);
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "Range table for loss for " << particleName << " is Retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
        if(fpi) {
          res = false;
	  G4cout << "Range table for loss for " << particleName << " from file <"
		 << filename << "> is not Retrieved"
		 << G4endl;
        }
      }

      filename = GetPhysicsTableFileName(part,directory,"DEDXnr",ascii);
      yes = theDEDXunRestrictedTable->ExistPhysicsTable(filename);
      if(yes) yes = G4PhysicsTableHelper::RetrievePhysicsTable(theDEDXunRestrictedTable,filename,ascii);
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "Non-restricted DEDX table for " << particleName << " is Retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
        if (1 < verboseLevel) {
          G4cout << "Non-restricted DEDX table for " << particleName << " from file <"
                 << filename << "> is not Retrieved"
                 << G4endl;
        }
      }

      filename = GetPhysicsTableFileName(part,directory,"PreciseRange",ascii);
      yes = thePreciseRangeTable->ExistPhysicsTable(filename);
      if(yes) yes = G4PhysicsTableHelper::RetrievePhysicsTable(thePreciseRangeTable,filename,ascii);
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "Precise Range table for " << particleName << " is Retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
	G4cout << "Precise Range table for loss for " << particleName << " does not exist"
	       << G4endl;
      }

      filename = GetPhysicsTableFileName(part,directory,"InverseRange",ascii);
      yes = theInverseRangeTable->ExistPhysicsTable(filename);
      if(yes)  yes = G4PhysicsTableHelper::RetrievePhysicsTable(theInverseRangeTable,filename,ascii);
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "InverseRange table for " << particleName << " is Retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
        if(fpi) {
          res = false;
          G4cout << "InverseRange table for " << particleName << " from file <"
                 << filename << "> is not Retrieved"
                 << G4endl;

        }
      }

      filename = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
      yes = theLambdaTable->ExistPhysicsTable(filename);
      if(yes) yes = G4PhysicsTableHelper::RetrievePhysicsTable(theLambdaTable,filename,ascii);
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "Lambda table for " << particleName << " is Retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
        if(fpi) {
          res = false;
          G4cout << "Lambda table for " << particleName << " from file <"
                 << filename << "> is not Retrieved"
                 << G4endl;
        }
      }

      filename = GetPhysicsTableFileName(part,directory,"SubLambda",ascii);
      yes = theSubLambdaTable->ExistPhysicsTable(filename);
      if(yes) yes = G4PhysicsTableHelper::RetrievePhysicsTable(theSubLambdaTable,filename,ascii);
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "SubLambda table for " << particleName << " is Retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
        if(nSCoffRegions) {
          res=false;
          G4cout << "SubLambda table for " << particleName << " from file <"
                 << filename << "> is not Retrieved"
                 << G4endl;
	}
      }
    }
  }

  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetLinearLossLimit(G4double val)
{
  linLossLimit = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetLossFluctuations(G4bool val)
{
  if(val && !lossFluctuationArePossible) return;
  lossFluctuationFlag = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetSubCutoff(G4bool)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetRandomStep(G4bool val)
{
  rndmStepFlag = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetMinSubRange(G4double val)
{
  minSubRange = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4VEnergyLossProcess::TablesAreBuilt() const
{
  return  tablesAreBuilt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4VEnergyLossProcess::NumberOfSubCutoffRegions() const
{
  return nSCoffRegions;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetDEDXBinning(G4int nbins)
{
  nDEDXBins = nbins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetDEDXBinningForPreciseRange(G4int nbins)
{
  nDEDXBinsForRange = nbins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetLambdaBinning(G4int nbins)
{
  nLambdaBins = nbins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::MinKinEnergy() const
{
  return minKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetMinKinEnergy(G4double e)
{
  minKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetMaxKinEnergy(G4double e)
{
  maxKinEnergy = e;
  if(e < maxKinEnergyForRange) maxKinEnergyForRange = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetMaxKinEnergyForPreciseRange(G4double e)
{
  maxKinEnergyForRange = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::MaxKinEnergy() const
{
  return maxKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::ActivateDeexcitation(G4bool, const G4Region*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetLambdaFactor(G4double val)
{
  if(val > 0.0 && val <= 1.0) lambdaFactor = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetIonisation(G4bool val)
{
  isIonisation = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4VEnergyLossProcess::IsIonisationProcess() const
{
  return isIonisation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

