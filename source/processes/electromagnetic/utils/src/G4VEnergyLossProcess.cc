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
// $Id: G4VEnergyLossProcess.cc,v 1.123 2008/01/11 19:55:29 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01-patch-01 $
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
// 25-07-05 Add extra protection PostStep for non-integral mode (V.Ivanchenko)
// 12-08-05 Integral=false; SetStepFunction(0.2, 0.1*mm) (mma)
// 18-08-05 Return back both AlongStep and PostStep from 7.0 (V.Ivanchenko)
// 02-09-05 Default StepFunction 0.2 1 mm + integral (V.Ivanchenko)
// 04-09-05 default lambdaFactor 0.8 (V.Ivanchenko)
// 05-10-05 protection against 0 energy loss added (L.Urban)
// 17-10-05 protection above has been removed (L.Urban)
// 06-01-06 reset currentCouple when StepFunction is changed (V.Ivanchenko)
// 10-01-06 PreciseRange -> CSDARange (V.Ivantchenko)
// 18-01-06 Clean up subcutoff including recalculation of presafety (VI)
// 20-01-06 Introduce G4EmTableType and reducing number of methods (VI)
// 22-03-06 Add control on warning printout AlongStep (VI)
// 23-03-06 Use isIonisation flag (V.Ivanchenko)
// 07-06-06 Do not reflect AlongStep in subcutoff regime (V.Ivanchenko)
// 14-01-07 add SetEmModel(index) and SetFluctModel() (mma)
// 16-01-07 add IonisationTable and IonisationSubTable (V.Ivanchenko)
// 16-02-07 set linLossLimit=1.e-6 (V.Ivanchenko)
// 13-03-07 use SafetyHelper instead of navigator (V.Ivanchenko)
// 10-04-07 use unique SafetyHelper (V.Ivanchenko)
// 12-04-07 Add verbosity at destruction (V.Ivanchenko)
// 25-04-07 move initialisation of safety helper to BuildPhysicsTable (VI)
// 27-10-07 Virtual functions moved to source (V.Ivanchenko)
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
#include "G4PhysicsLogVector.hh"
#include "G4VParticleChange.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4ProcessManager.hh"
#include "G4UnitsTable.hh"
#include "G4GenericIon.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4SafetyHelper.hh"
#include "G4TransportationManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossProcess::G4VEnergyLossProcess(const G4String& name, 
  G4ProcessType type): G4VContinuousDiscreteProcess(name, type),
  nSCoffRegions(0),
  idxSCoffRegions(0),
  nProcesses(0),
  theDEDXTable(0),
  theDEDXSubTable(0),
  theDEDXunRestrictedTable(0),
  theIonisationTable(0),
  theIonisationSubTable(0),
  theRangeTableForLoss(0),
  theCSDARangeTable(0),
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
  nBins(120),
  nBinsCSDA(70),
  nWarnings(0),
  linLossLimit(0.05),
  minSubRange(0.1),
  lambdaFactor(0.8),
  mfpKinEnergy(0.0),
  lossFluctuationFlag(true),
  rndmStepFlag(false),
  tablesAreBuilt(false),
  integral(true),
  isIonisation(true),
  useSubCutoff(false)
{
  SetVerboseLevel(1);

  // Size of tables
  lowestKinEnergy  = 1.*eV;
  minKinEnergy     = 0.1*keV;
  maxKinEnergy     = 100.0*TeV;
  maxKinEnergyCSDA = 1.0*GeV;

  // default dRoverRange and finalRange
  SetStepFunction(0.2, 1.0*mm);

  theElectron = G4Electron::Electron();
  thePositron = G4Positron::Positron();

  // run time objects
  pParticleChange = &fParticleChange;
  modelManager = new G4EmModelManager();
  safetyHelper = G4TransportationManager::GetTransportationManager()
    ->GetSafetyHelper();
  aGPILSelection = CandidateForSelection;

  // initialise model
  (G4LossTableManager::Instance())->Register(this);
  fluctModel = 0;

  scoffRegions.clear();
  scProcesses.clear();
  scTracks.reserve(5);
  secParticles.reserve(5);

  // Data for stragling of ranges from ICRU'37 report
  const G4int nrbins = 7;
  vstrag = new G4PhysicsLogVector(keV, GeV, nrbins);
  G4double s[nrbins] = {-0.2, -0.85, -1.3, -1.578, -1.76, -1.85, -1.9};
  for(G4int i=0; i<nrbins; i++) {vstrag->PutValue(i, s[i]);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossProcess::~G4VEnergyLossProcess()
{
  if(1 < verboseLevel) 
    G4cout << "G4VEnergyLossProcess destruct " << GetProcessName() 
	   << G4endl;
  delete vstrag;
  Clear();

  if ( !baseParticle ) {
    if(theDEDXTable && theRangeTableForLoss) {
      if(theIonisationTable == theDEDXTable) theIonisationTable = 0;
      theDEDXTable->clearAndDestroy();
      delete theDEDXTable;
      if(theDEDXSubTable) {
	if(theIonisationSubTable == theDEDXSubTable) 
	  theIonisationSubTable = 0;
	theDEDXSubTable->clearAndDestroy();
        delete theDEDXSubTable;
      }
    }
    if(theIonisationTable) {
      theIonisationTable->clearAndDestroy(); 
      delete theIonisationTable;
    }
    if(theIonisationSubTable) {
      theIonisationSubTable->clearAndDestroy(); 
      delete theIonisationSubTable;
    }
    if(theDEDXunRestrictedTable && theCSDARangeTable) {
       theDEDXunRestrictedTable->clearAndDestroy();
       delete theDEDXunRestrictedTable;
    }
    if(theCSDARangeTable) {
      theCSDARangeTable->clearAndDestroy();
      delete theCSDARangeTable;
    }
    if(theRangeTableForLoss) {
      theRangeTableForLoss->clearAndDestroy();
      delete theRangeTableForLoss;
    }
    if(theInverseRangeTable) {
      theInverseRangeTable->clearAndDestroy();
      delete theInverseRangeTable;
    }
    if(theLambdaTable) {
      theLambdaTable->clearAndDestroy();
      delete theLambdaTable;
    }
    if(theSubLambdaTable) {
      theSubLambdaTable->clearAndDestroy();
      delete theSubLambdaTable;
    }
  }

  delete modelManager;
  (G4LossTableManager::Instance())->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::Clear()
{
  if(1 < verboseLevel) 
    G4cout << "G4VEnergyLossProcess::Clear() for " << GetProcessName() 
    << G4endl;

  delete [] theDEDXAtMaxEnergy;
  delete [] theRangeAtMaxEnergy;
  delete [] theEnergyOfCrossSectionMax;
  delete [] theCrossSectionMax;
  delete [] idxSCoffRegions;

  theDEDXAtMaxEnergy = 0;
  theRangeAtMaxEnergy = 0;
  theEnergyOfCrossSectionMax = 0;
  theCrossSectionMax = 0;
  tablesAreBuilt = false;

  scTracks.clear();
  scProcesses.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::PreparePhysicsTable(
     const G4ParticleDefinition& part)
{

  // Are particle defined?
  if( !particle ) {
    if(part.GetParticleType() == "nucleus" && 
       part.GetParticleSubType() == "generic") 
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

  if(&part != particle) {
    if(part.GetParticleType() == "nucleus") lManager->RegisterIon(&part, this);
    else                          lManager->RegisterExtraParticle(&part, this);
    return;
  }

  Clear();

  currentCouple = 0;
  preStepLambda = 0.0;
  mfpKinEnergy  = DBL_MAX;
  fRange        = DBL_MAX;
  preStepKinEnergy = 0.0;

  // Base particle and set of models can be defined here
  InitialiseEnergyLossProcess(particle, baseParticle);

  // Tables preparation
  if (!baseParticle) {
    
    theDEDXTable = G4PhysicsTableHelper::PreparePhysicsTable(theDEDXTable);
    if (lManager->BuildCSDARange()) {
      theDEDXunRestrictedTable = 
	G4PhysicsTableHelper::PreparePhysicsTable(theDEDXunRestrictedTable);
      theCSDARangeTable = 
	G4PhysicsTableHelper::PreparePhysicsTable(theCSDARangeTable);
    }

    theRangeTableForLoss = 
      G4PhysicsTableHelper::PreparePhysicsTable(theRangeTableForLoss);
    theInverseRangeTable = 
      G4PhysicsTableHelper::PreparePhysicsTable(theInverseRangeTable);
  
    theLambdaTable = G4PhysicsTableHelper::PreparePhysicsTable(theLambdaTable);
    if (nSCoffRegions) {
      theDEDXSubTable = 
	G4PhysicsTableHelper::PreparePhysicsTable(theDEDXSubTable);
      theSubLambdaTable = 
	G4PhysicsTableHelper::PreparePhysicsTable(theSubLambdaTable);
    }
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
    if(chargeSqRatio > 0.0) reduceFactor = 1.0/(chargeSqRatio*massRatio);
  }

  theCuts = modelManager->Initialise(particle, secondaryParticle, 
				     minSubRange, verboseLevel);

  // Sub Cutoff Regime
  scProcesses.clear();
  nProcesses = 0;
  
  if (nSCoffRegions>0) {
    theSubCuts = modelManager->SubCutoff();

    const G4ProductionCutsTable* theCoupleTable=
          G4ProductionCutsTable::GetProductionCutsTable();
    size_t numOfCouples = theCoupleTable->GetTableSize();
    idxSCoffRegions = new G4int[numOfCouples];
  
    for (size_t j=0; j<numOfCouples; j++) {

      const G4MaterialCutsCouple* couple = 
	theCoupleTable->GetMaterialCutsCouple(j);
      const G4ProductionCuts* pcuts = couple->GetProductionCuts();
      G4int reg = 0;
      for(G4int i=0; i<nSCoffRegions; i++) {
        if( pcuts == scoffRegions[i]->GetProductionCuts()) reg = 1;
      }
      idxSCoffRegions[j] = reg;
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
    G4cout << "### G4VEnergyLossProcess::BuildPhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << "; local: " << particle->GetParticleName();
    if(baseParticle) G4cout << "; base: " << baseParticle->GetParticleName();
    G4cout << G4endl;
  }

  if(!tablesAreBuilt && &part == particle)
    G4LossTableManager::Instance()->BuildPhysicsTable(particle, this);

  if(0 < verboseLevel && (&part == particle) && !baseParticle) {
    PrintInfoDefinition();
    safetyHelper->InitialiseHelper();
  }

  if(1 < verboseLevel) {
    G4cout << "### G4VEnergyLossProcess::BuildPhysicsTable() done for "
           << GetProcessName()
           << " and particle " << part.GetParticleName();
    if(isIonisation) G4cout << "  isIonisation  flag = 1";
    G4cout << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::ActivateSubCutoff(G4bool val, const G4Region* r)
{
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  if(val) {
    useSubCutoff = true;
    if (!r) r = regionStore->GetRegion("DefaultRegionForTheWorld", false);
    if (nSCoffRegions) {
      for (G4int i=0; i<nSCoffRegions; i++) {
	if (r == scoffRegions[i]) return;
      }
    }
    scoffRegions.push_back(r);
    nSCoffRegions++;
  } else {
    useSubCutoff = false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossProcess::BuildDEDXTable(G4EmTableType tType)
{
  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTable() of type " << tType
	   << " for " << GetProcessName()
           << " and particle " << particle->GetParticleName()
           << G4endl;
  }
  G4PhysicsTable* table = 0;
  G4double emin = minKinEnergy;
  G4double emax = maxKinEnergy;
  G4int bin = nBins;

  if(fTotal == tType) {
    emax  = maxKinEnergyCSDA;
    bin   = nBinsCSDA;
    table = theDEDXunRestrictedTable;
  } else if(fRestricted == tType) {
    table = theDEDXTable;
    if(theIonisationTable) 
      table = G4PhysicsTableHelper::PreparePhysicsTable(theIonisationTable); 
  } else if(fSubRestricted == tType) {    
    table = theDEDXSubTable;
    if(theIonisationSubTable) 
      table = G4PhysicsTableHelper::PreparePhysicsTable(theIonisationSubTable); 
  } else {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTable WARNING: wrong type "
	   << tType << G4endl;
  }

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  if(1 < verboseLevel) {
    G4cout << numOfCouples << " materials"
           << " minKinEnergy= " << minKinEnergy
           << " maxKinEnergy= " << maxKinEnergy
           << " EmTableType= " << tType
           << " table= " << table
           << G4endl;
  }
  if(!table) return table;

  for(size_t i=0; i<numOfCouples; i++) {

    if(1 < verboseLevel) 
      G4cout << "G4VEnergyLossProcess::BuildDEDXVector flag=  " 
	     << table->GetFlag(i) << G4endl;

    if (table->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = 
	theCoupleTable->GetMaterialCutsCouple(i);
      G4PhysicsVector* aVector = new G4PhysicsLogVector(emin, emax, bin);
      modelManager->FillDEDXVector(aVector, couple, tType);

      // Insert vector for this material into the table
      G4PhysicsTableHelper::SetPhysicsVector(table, i, aVector);
    }
  }

  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTable(): table is built for "
           << particle->GetParticleName()
           << " and process " << GetProcessName()
           << G4endl;
    //    if(2 < verboseLevel) G4cout << (*table) << G4endl;
  }

  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossProcess::BuildLambdaTable(G4EmTableType tType)
{
  G4PhysicsTable* table = 0;

  if(fRestricted == tType) {
    table = theLambdaTable;
  } else if(fSubRestricted == tType) {    
    table = theSubLambdaTable;
  } else {
    G4cout << "G4VEnergyLossProcess::BuildLambdaTable WARNING: wrong type "
	   << tType << G4endl;
  }

  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildLambdaTable() of type "
	   << tType << " for process "
           << GetProcessName() << " and particle "
           << particle->GetParticleName()
           << " EmTableType= " << tType
           << " table= " << table
           << G4endl;
  }
  if(!table) return table;

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  for(size_t i=0; i<numOfCouples; i++) {

    if (table->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = 
	theCoupleTable->GetMaterialCutsCouple(i);
      G4double cut = (*theCuts)[i];
      if(fSubRestricted == tType) cut = (*theSubCuts)[i]; 
      G4PhysicsVector* aVector = LambdaPhysicsVector(couple, cut);
      modelManager->FillLambdaVector(aVector, couple, true, tType);

      // Insert vector for this material into the table
      G4PhysicsTableHelper::SetPhysicsVector(table, i, aVector);
    }
  }

  if(1 < verboseLevel) {
    G4cout << "Lambda table is built for "
           << particle->GetParticleName()
           << G4endl;
  }

  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::GetContinuousStepLimit(
		const G4Track&,
                G4double, G4double, G4double&)
{
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::AlongStepGetPhysicalInteractionLength(
                             const G4Track&,
                             G4double,
                             G4double  currentMinStep,
                             G4double&,
                             G4GPILSelection* selection)
{
  G4double x = DBL_MAX;
  *selection = aGPILSelection;
  if(isIonisation) {
    fRange = GetScaledRangeForScaledEnergy(preStepScaledEnergy)*reduceFactor;

    x = fRange;
    G4double y = x*dRoverRange;

    if(x > finalRange && y < currentMinStep) { 
      x = y + finalRange*(1.0 - dRoverRange)*(2.0 - finalRange/fRange);
    } else if (rndmStepFlag) x = SampleRange();
    //    G4cout<<GetProcessName()<<": e= "<<preStepKinEnergy
    //	  <<" range= "<<fRange <<" cMinSt="<<currentMinStep
    //	  <<" safety= " << safety<< " limit= " << x <<G4endl;
  }
  //  G4cout<<GetProcessName()<<": e= "<<preStepKinEnergy
  //  <<" stepLimit= "<<x<<G4endl;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::GetMeanFreePath(
                             const G4Track& track,
                             G4double,
                             G4ForceCondition* condition)

{
  *condition = NotForced;
  return MeanFreePath(track);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  G4double x = DBL_MAX;
  if(previousStepSize <= DBL_MIN) theNumberOfInteractionLengthLeft = -1.0;
  InitialiseStep(track);

  if(preStepScaledEnergy < mfpKinEnergy) {
    if (integral) ComputeLambdaForScaledEnergy(preStepScaledEnergy);
    else  preStepLambda = GetLambdaForScaledEnergy(preStepScaledEnergy);
    if(preStepLambda <= DBL_MIN) mfpKinEnergy = 0.0;
  }

  // non-zero cross section
  if(preStepLambda > DBL_MIN) { 
    if (theNumberOfInteractionLengthLeft < 0.0) {
      // beggining of tracking (or just after DoIt of this process)
      ResetNumberOfInteractionLengthLeft();
    } else if(currentInteractionLength < DBL_MAX) {
      // subtract NumberOfInteractionLengthLeft
      SubtractNumberOfInteractionLengthLeft(previousStepSize);
      if(theNumberOfInteractionLengthLeft < 0.)
	theNumberOfInteractionLengthLeft = perMillion;
    }

    // get mean free path and step limit
    currentInteractionLength = 1.0/preStepLambda;
    x = theNumberOfInteractionLengthLeft * currentInteractionLength;

#ifdef G4VERBOSE
    if (verboseLevel>2){
      G4cout << "G4VEnergyLossProcess::PostStepGetPhysicalInteractionLength ";
      G4cout << "[ " << GetProcessName() << "]" << G4endl; 
      G4cout << " for " << particle->GetParticleName() 
             << " in Material  " <<  currentMaterial->GetName()
	     << " Ekin(MeV)= " << preStepKinEnergy/MeV 
	     <<G4endl;
      G4cout << "MeanFreePath = " << currentInteractionLength/cm << "[cm]" 
	     << "InteractionLength= " << x/cm <<"[cm] " <<G4endl;
    }
#endif

    // zero cross section case
  } else {
    if(theNumberOfInteractionLengthLeft > DBL_MIN && 
       currentInteractionLength < DBL_MAX) {
      // subtract NumberOfInteractionLengthLeft
      SubtractNumberOfInteractionLengthLeft(previousStepSize);
      if(theNumberOfInteractionLengthLeft < 0.)
	theNumberOfInteractionLengthLeft = perMillion;
    }
    currentInteractionLength = DBL_MAX;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEnergyLossProcess::AlongStepDoIt(const G4Track& track,
                                                       const G4Step& step)
{
  fParticleChange.InitializeForAlongStep(track);
  // The process has range table - calculate energy loss
  if(!isIonisation) return &fParticleChange;

  // Get the actual (true) Step length
  G4double length = step.GetStepLength();
  if(length <= DBL_MIN) return &fParticleChange;
  G4double eloss  = 0.0;

  /*
  if(-1 < verboseLevel) {
    const G4ParticleDefinition* d = track.GetDefinition();
    G4cout << "AlongStepDoIt for "
           << GetProcessName() << " and particle "
           << d->GetParticleName()
           << "  eScaled(MeV)= " << preStepScaledEnergy/MeV
           << "  range(mm)= " << fRange/mm
           << "  s(mm)= " << length/mm
           << "  q^2= " << chargeSqRatio
           << " md= " << d->GetPDGMass()
           << "  status= " << track.GetTrackStatus()
           << G4endl;
  }
  */

  // stopping
  if (length >= fRange) {
    fParticleChange.SetProposedKineticEnergy(0.0);
    fParticleChange.ProposeLocalEnergyDeposit(preStepKinEnergy);
    return &fParticleChange;
  }

  // Short step
  eloss = GetDEDXForScaledEnergy(preStepScaledEnergy)*length;

  // Long step
  //} else {
  if(eloss > preStepKinEnergy*linLossLimit) {

    G4double x = 
      GetScaledRangeForScaledEnergy(preStepScaledEnergy) - length/reduceFactor;
    eloss = preStepKinEnergy - ScaledKinEnergyForLoss(x)/massRatio;
    /*
    if(-1 < verboseLevel) 
      G4cout << "Long STEP: rPre(mm)= " 
             << GetScaledRangeForScaledEnergy(preStepScaledEnergy)/mm
             << " rPost(mm)= " << x/mm
             << " ePre(MeV)= " << preStepScaledEnergy/MeV
             << " eloss(MeV)= " << eloss/MeV
             << " eloss0(MeV)= "
             << GetDEDXForScaledEnergy(preStepScaledEnergy)*length/MeV
             << G4endl;
    */
  }

  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  G4VEmModel* currentModel = SelectModel(preStepScaledEnergy);
  /*    
  G4double eloss0 = eloss;
  if(-1 < verboseLevel ) {
    G4cout << "Before fluct: eloss(MeV)= " << eloss/MeV
           << " e-eloss= " << preStepKinEnergy-eloss
           << " step(mm)= " << length/mm
           << " range(mm)= " << fRange/mm
           << " fluct= " << lossFluctuationFlag
           << G4endl;
  }
  */

  G4double cut  = (*theCuts)[currentMaterialIndex];
  G4double esec = 0.0;
  G4double esecdep = 0.0;

  // SubCutOff 
  if(useSubCutoff) {
    if(idxSCoffRegions[currentMaterialIndex]) {

      G4double currentMinSafety = 0.0; 
      G4StepPoint* prePoint  = step.GetPreStepPoint();
      G4StepPoint* postPoint = step.GetPostStepPoint();
      G4double preSafety  = prePoint->GetSafety();
      G4double postSafety = preSafety - length; 
      G4double rcut = currentCouple->GetProductionCuts()->GetProductionCut(1);

      // recompute safety
      if(prePoint->GetStepStatus() != fGeomBoundary &&
	 postPoint->GetStepStatus() != fGeomBoundary) {

	/*
	//      G4bool yes = (track.GetTrackID() == 5512);
        G4bool yes = false;
	if(yes)
	  G4cout << "G4VEnergyLoss: presafety= " << preSafety
		 << " rcut= " << rcut << "  length= " << length 
		 << " dir " << track.GetMomentumDirection()
		 << G4endl;
	*/

	if(preSafety < rcut) 
	  preSafety = safetyHelper->ComputeSafety(prePoint->GetPosition());

	//if(yes) {
	//  G4cout << "G4VEnergyLoss: newsafety= " << preSafety << G4endl;
	  //	   if(preSafety==0.0 && track.GetTrackID() == 5512 ) exit(1);
	//}
	if(postSafety < rcut) 
	  postSafety = safetyHelper->ComputeSafety(postPoint->GetPosition());
	/*	
	  if(-1 < verboseLevel) 
	  G4cout << "Subcutoff: presafety(mm)= " << preSafety/mm
	         << " postsafety(mm)= " << postSafety/mm
	         << " rcut(mm)= " << rcut/mm 
	         << G4endl;
	*/
	currentMinSafety = std::min(preSafety,postSafety); 
      }

      // Decide to start subcut sampling
      if(currentMinSafety < rcut) {

        cut = (*theSubCuts)[currentMaterialIndex];
 	eloss -= GetSubDEDXForScaledEnergy(preStepScaledEnergy)*length;
	scTracks.clear();
	SampleSubCutSecondaries(scTracks, step, 
				currentModel,currentMaterialIndex, 
				esecdep);
	/*
	if(nProcesses > 0) {
	  for(G4int i=0; i<nProcesses; i++) {
	    (scProcesses[i])->SampleSubCutSecondaries(
		scTracks, step, (scProcesses[i])->
		SelectModelForMaterial(preStepKinEnergy, currentMaterialIndex),
		currentMaterialIndex,esecdep);
	  }
	} 
	*/   
	G4int n = scTracks.size();
	if(n>0) {
	  G4ThreeVector mom = dynParticle->GetMomentum();
	  fParticleChange.SetNumberOfSecondaries(n);
	  for(G4int i=0; i<n; i++) {
	    G4Track* t = scTracks[i];
	    G4double e = t->GetKineticEnergy();
	    if (t->GetDefinition() == thePositron) e += 2.0*electron_mass_c2;
	    esec += e;
	    pParticleChange->AddSecondary(t);
	    //  mom -= t->GetMomentum();
	  }      
	  //	    fParticleChange.SetProposedMomentum(mom);            
	}
      }
    }
  }

  // Corrections, which cannot be tabulated
  CorrectionsAlongStep(currentCouple, dynParticle, eloss, length);

  // Sample fluctuations
  if (lossFluctuationFlag) {
    G4VEmFluctuationModel* fluc = currentModel->GetModelOfFluctuations();
    if(fluc && 
      (eloss + esec + esecdep + lowestKinEnergy) < preStepKinEnergy) {

      G4double tmax = 
	std::min(currentModel->MaxSecondaryKinEnergy(dynParticle),cut);
      eloss = fluc->SampleFluctuations(currentMaterial,dynParticle,
				       tmax,length,eloss);
      /*           
      if(-1 < verboseLevel) 
      G4cout << "After fluct: eloss(MeV)= " << eloss/MeV
             << " fluc= " << (eloss-eloss0)/MeV
             << " currentChargeSquare= " << chargeSquare
             << " massRatio= " << massRatio
             << " tmax= " << tmax
             << G4endl;
      */
    }
  }
  // add low-energy subcutoff particles
  eloss += esecdep;
  if(eloss < 0.0) eloss = 0.0;

  // Energy balanse
  G4double finalT = preStepKinEnergy - eloss - esec;
  if (finalT <= lowestKinEnergy) {
    eloss  = preStepKinEnergy - esec;
    finalT = 0.0;
  }

  fParticleChange.SetProposedKineticEnergy(finalT);
  fParticleChange.ProposeLocalEnergyDeposit(eloss);

  /*
  if(-1 < verboseLevel) {
    G4cout << "Final value eloss(MeV)= " << eloss/MeV
           << " preStepKinEnergy= " << preStepKinEnergy
           << " postStepKinEnergy= " << finalT
           << " lossFlag= " << lossFluctuationFlag
           << "  status= " << track.GetTrackStatus()
           << G4endl;
  }
  */

  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SampleSubCutSecondaries(
       std::vector<G4Track*>& tracks, 
       const G4Step& step, 
       G4VEmModel* model,
       G4int idx,
       G4double& extraEdep) 
{
  // Fast check weather subcutoff can work
  G4double subcut = (*theSubCuts)[idx];
  G4double cut = (*theCuts)[idx];
  if(cut <= subcut) return;

  const G4Track* track = step.GetTrack();
  const G4DynamicParticle* dp = track->GetDynamicParticle();
  G4bool b;
  G4double cross = 
    chargeSqRatio*(((*theSubLambdaTable)[idx])->GetValue(dp->GetKineticEnergy(),b));
  G4double length = step.GetStepLength();

  // negligible probability to get any interaction
  if(length*cross < perMillion) return;
  /*    
  if(-1 < verboseLevel) 
    G4cout << "<<< Subcutoff for " << GetProcessName()
	   << " cross(1/mm)= " << cross*mm << ">>>"
	   << " e(MeV)= " << preStepScaledEnergy
	   << " matIdx= " << currentMaterialIndex
	   << G4endl;
  */

  // Sample subcutoff secondaries
  G4StepPoint* preStepPoint = step.GetPreStepPoint();
  G4ThreeVector prepoint = preStepPoint->GetPosition();
  G4ThreeVector dr = step.GetPostStepPoint()->GetPosition() - prepoint;
  G4double pretime = preStepPoint->GetGlobalTime();
  //  G4double dt = length/preStepPoint->GetVelocity();
  G4double fragment = 0.0;

  do {
    G4double del = -std::log(G4UniformRand())/cross;
    fragment += del/length;
    if (fragment > 1.0) break;

    // sample secondaries
    secParticles.clear();
    model->SampleSecondaries(&secParticles,track->GetMaterialCutsCouple(),
			     dp,subcut,cut);

    // position of subcutoff particles
    G4ThreeVector r = prepoint + fragment*dr;
    std::vector<G4DynamicParticle*>::iterator it;
    for(it=secParticles.begin(); it!=secParticles.end(); it++) {

      G4bool addSec = true;
      // do not track very low-energy delta-electrons
      if(theSecondaryRangeTable && (*it)->GetDefinition() == theElectron) {
	G4bool b;
	G4double ekin = (*it)->GetKineticEnergy();
	G4double rg = ((*theSecondaryRangeTable)[idx]->GetValue(ekin, b));
	//          if(rg < currentMinSafety) {
	if(rg < safetyHelper->ComputeSafety(r)) {
	  extraEdep += ekin;
	  delete (*it);
	  addSec = false;
	}
      }
      if(addSec) {
	//	G4Track* t = new G4Track((*it), pretime + fragment*dt, r);
	G4Track* t = new G4Track((*it), pretime, r);
	t->SetTouchableHandle(track->GetTouchableHandle());
	tracks.push_back(t);

	/*
	  if(-1 < verboseLevel) 
	  G4cout << "New track " << p->GetDefinition()->GetParticleName()
	  << " e(keV)= " << p->GetKineticEnergy()/keV
	  << " fragment= " << fragment
	  << G4endl;
	*/
      }
    }
  } while (fragment <= 1.0);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEnergyLossProcess::PostStepDoIt(const G4Track& track,
                                                      const G4Step&)
{
  fParticleChange.InitializeForPostStep(track);
  G4double finalT = track.GetKineticEnergy();
  if(finalT <= lowestKinEnergy) return &fParticleChange;

  G4double postStepScaledEnergy = finalT*massRatio;
  /*
  if(-1 < verboseLevel) {
    G4cout << GetProcessName()
           << "::PostStepDoIt: E(MeV)= " << finalT/MeV
	   << G4endl;
  }
  */
  // Integral approach
  if (integral) {
    G4double lx = GetLambdaForScaledEnergy(postStepScaledEnergy);
    /*
    if(preStepLambda<lx && 1 < verboseLevel && nWarnings<200) {
      G4cout << "WARNING: for " << particle->GetParticleName()
             << " and " << GetProcessName()
             << " E(MeV)= " << finalT/MeV
             << " preLambda= " << preStepLambda 
	     << " < " << lx << " (postLambda) "
	     << G4endl;
      nWarnings++;
    }
    */
    if(preStepLambda*G4UniformRand() > lx) {
      ClearNumberOfInteractionLengthLeft();
      return &fParticleChange;
    }
  }

  G4VEmModel* currentModel = SelectModel(postStepScaledEnergy);
  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  G4double tcut = (*theCuts)[currentMaterialIndex];

  // sample secondaries
  secParticles.clear();
  currentModel->SampleSecondaries(&secParticles, currentCouple, dynParticle, tcut);

  // save secondaries
  G4int num = secParticles.size();
  if(num > 0) {
    fParticleChange.SetNumberOfSecondaries(num);
    for (G4int i=0; i<num; i++) {
      fParticleChange.AddSecondary(secParticles[i]);
    }
  }

  /*
  if(-1 < verboseLevel) {
    G4cout << "::PostStepDoIt: Sample secondary; Efin= " 
    << fParticleChange.GetProposedKineticEnergy()/MeV
           << " MeV; model= (" << currentModel->LowEnergyLimit()
           << ", " <<  currentModel->HighEnergyLimit() << ")"
           << "  preStepLambda= " << preStepLambda
           << "  dir= " << track.GetMomentumDirection()
           << "  status= " << track.GetTrackStatus()
           << G4endl;
  }
  */
  ClearNumberOfInteractionLengthLeft();
  return &fParticleChange;
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
           << " in " << nBins << " bins." << G4endl
           << "      Lambda tables from threshold to "
           << G4BestUnit(maxKinEnergy,"Energy")
           << " in " << nBins << " bins."
           << G4endl;
    PrintInfo();
    if(theRangeTableForLoss && isIonisation) 
      G4cout << "      Step function: finalRange(mm)= " << finalRange/mm
             << ", dRoverRange= " << dRoverRange
             << ", integral: " << integral
             << ", fluct: " << lossFluctuationFlag
             << G4endl;
    
    if(theCSDARangeTable && isIonisation) 
      G4cout << "      CSDA range table up"
             << " to " << G4BestUnit(maxKinEnergyCSDA,"Energy")
             << " in " << nBinsCSDA << " bins." << G4endl;
    
    if(nSCoffRegions>0) 
      G4cout << "      Subcutoff sampling in " << nSCoffRegions 
	     << " regions" << G4endl;

    if(2 < verboseLevel) {
      G4cout << "DEDXTable address= " << theDEDXTable << G4endl;
      if(theDEDXTable && isIonisation) G4cout << (*theDEDXTable) << G4endl;
      G4cout << "non restricted DEDXTable address= " 
	     << theDEDXunRestrictedTable << G4endl;
      if(theDEDXunRestrictedTable && isIonisation) 
           G4cout << (*theDEDXunRestrictedTable) << G4endl;
      if(theDEDXSubTable && isIonisation) G4cout << (*theDEDXSubTable) 
						 << G4endl;
      G4cout << "CSDARangeTable address= " << theCSDARangeTable 
	     << G4endl;
      if(theCSDARangeTable && isIonisation) G4cout << (*theCSDARangeTable) 
            << G4endl;
      G4cout << "RangeTableForLoss address= " << theRangeTableForLoss 
	     << G4endl;
      if(theRangeTableForLoss && isIonisation) 
             G4cout << (*theRangeTableForLoss) << G4endl;
      G4cout << "InverseRangeTable address= " << theInverseRangeTable 
	     << G4endl;
      if(theInverseRangeTable && isIonisation) 
             G4cout << (*theInverseRangeTable) << G4endl;
      G4cout << "LambdaTable address= " << theLambdaTable << G4endl;
      if(theLambdaTable && isIonisation) G4cout << (*theLambdaTable) << G4endl;
      G4cout << "SubLambdaTable address= " << theSubLambdaTable << G4endl;
      if(theSubLambdaTable && isIonisation) G4cout << (*theSubLambdaTable) 
             << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetDEDXTable(G4PhysicsTable* p, G4EmTableType tType)
{
  if(fTotal == tType && theDEDXunRestrictedTable != p) {
    if(theDEDXunRestrictedTable) theDEDXunRestrictedTable->clearAndDestroy();
    theDEDXunRestrictedTable = p;
    if(p) {
      size_t n = p->length();
      G4PhysicsVector* pv = (*p)[0];
      G4double emax = maxKinEnergyCSDA;
      G4bool b;
      theDEDXAtMaxEnergy = new G4double [n];

      for (size_t i=0; i<n; i++) {
	pv = (*p)[i];
	G4double dedx = pv->GetValue(emax, b);
	theDEDXAtMaxEnergy[i] = dedx;
	//G4cout << "i= " << i << " emax(MeV)= " << emax/MeV<< " dedx= " 
	//<< dedx << G4endl;
      }
    }

  } else if(fRestricted == tType) {
    theDEDXTable = p;
  } else if(fSubRestricted == tType) {    
    theDEDXSubTable = p;
  } else if(fIonisation == tType && theIonisationTable != p) {    
    if(theIonisationTable) theIonisationTable->clearAndDestroy();
    theIonisationTable = p;
  } else if(fSubIonisation == tType && theIonisationSubTable != p) {    
    if(theIonisationSubTable) theIonisationSubTable->clearAndDestroy();
    theIonisationSubTable = p;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetCSDARangeTable(G4PhysicsTable* p)
{
  if(theCSDARangeTable != p) theCSDARangeTable = p;

  if(p) {
    size_t n = p->length();
    G4PhysicsVector* pv = (*p)[0];
    G4double emax = maxKinEnergyCSDA;
    G4bool b;
    theRangeAtMaxEnergy = new G4double [n];

    for (size_t i=0; i<n; i++) {
      pv = (*p)[i];
      G4double r2 = pv->GetValue(emax, b);
      theRangeAtMaxEnergy[i] = r2;
      //G4cout << "i= " << i << " e2(MeV)= " << emax/MeV << " r2= " 
      //<< r2<< G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetRangeTableForLoss(G4PhysicsTable* p)
{
  if(theRangeTableForLoss != p) {
    theRangeTableForLoss = p;
    if(1 < verboseLevel) {
      G4cout << "### Set Range table " << p 
	     << " for " << particle->GetParticleName()
             << " and process " << GetProcessName() << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetSecondaryRangeTable(G4PhysicsTable* p)
{
  if(theSecondaryRangeTable != p) {
    theSecondaryRangeTable = p;
    if(1 < verboseLevel) {
      G4cout << "### Set SecondaryRange table " << p 
	     << " for " << particle->GetParticleName()
             << " and process " << GetProcessName() << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetInverseRangeTable(G4PhysicsTable* p)
{
  if(theInverseRangeTable != p) {
    theInverseRangeTable = p;
    if(1 < verboseLevel) {
      G4cout << "### Set InverseRange table " << p 
	     << " for " << particle->GetParticleName()
             << " and process " << GetProcessName() << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetLambdaTable(G4PhysicsTable* p)
{
  if(1 < verboseLevel) {
    G4cout << "### Set Lambda table " << p 
	   << " for " << particle->GetParticleName()
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
  if(theSubLambdaTable != p) {
    theSubLambdaTable = p;
    if(1 < verboseLevel) {
      G4cout << "### Set SebLambda table " << p 
	     << " for " << particle->GetParticleName()
             << " and process " << GetProcessName() << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEnergyLossProcess::LambdaPhysicsVector(
                 const G4MaterialCutsCouple* couple, G4double cut)
{
  //  G4double cut  = (*theCuts)[couple->GetIndex()];
  //  G4int nbins = nLambdaBins;
  G4double tmin = 
    std::max(MinPrimaryEnergy(particle, couple->GetMaterial(), cut),
	     minKinEnergy);
  if(tmin >= maxKinEnergy) tmin = 0.5*maxKinEnergy;
  G4PhysicsVector* v = new G4PhysicsLogVector(tmin, maxKinEnergy, nBins);
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::MicroscopicCrossSection(
         G4double kineticEnergy, const G4MaterialCutsCouple* couple)
{
  // Cross section per atom is calculated
  DefineMaterial(couple);
  G4double cross = 0.0;
  G4bool b;
  if(theLambdaTable) 
    cross = 
      ((*theLambdaTable)[currentMaterialIndex])->GetValue(kineticEnergy, b)/
      currentMaterial->GetTotNbOfAtomsPerVolume();

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4VEnergyLossProcess::StorePhysicsTable(
       const G4ParticleDefinition* part, const G4String& directory, 
       G4bool ascii)
{
  G4bool res = true;
  if ( baseParticle || part != particle ) return res;

  if ( theDEDXTable ) {
    const G4String name = GetPhysicsTableFileName(part,directory,"DEDX",ascii);
    if( !theDEDXTable->StorePhysicsTable(name,ascii)) res = false;
  }

  if ( theDEDXunRestrictedTable ) {
    const G4String name = 
      GetPhysicsTableFileName(part,directory,"DEDXnr",ascii);
    if( !theDEDXTable->StorePhysicsTable(name,ascii)) res = false;
  }

  if ( theDEDXSubTable ) {
    const G4String name = 
      GetPhysicsTableFileName(part,directory,"SubDEDX",ascii);
    if( !theDEDXSubTable->StorePhysicsTable(name,ascii)) res = false;
  }

  if ( theIonisationTable ) {
    const G4String name = 
      GetPhysicsTableFileName(part,directory,"Ionisation",ascii);
    if( !theIonisationTable->StorePhysicsTable(name,ascii)) res = false;
  }

  if ( theIonisationSubTable ) {
    const G4String name = 
      GetPhysicsTableFileName(part,directory,"SubIonisation",ascii);
    if( !theIonisationSubTable->StorePhysicsTable(name,ascii)) res = false;
  }

  if ( theCSDARangeTable && isIonisation ) {
    const G4String name = 
      GetPhysicsTableFileName(part,directory,"CSDARange",ascii);
    if( !theCSDARangeTable->StorePhysicsTable(name,ascii)) res = false;
  }

  if ( theRangeTableForLoss && isIonisation ) {
    const G4String name = 
      GetPhysicsTableFileName(part,directory,"Range",ascii);
    if( !theRangeTableForLoss->StorePhysicsTable(name,ascii)) res = false;
  }

  if ( theInverseRangeTable && isIonisation ) {
    const G4String name = 
      GetPhysicsTableFileName(part,directory,"InverseRange",ascii);
    if( !theInverseRangeTable->StorePhysicsTable(name,ascii)) res = false;
  }

  if ( theLambdaTable  && isIonisation) {
    const G4String name = 
      GetPhysicsTableFileName(part,directory,"Lambda",ascii);
    if( !theLambdaTable->StorePhysicsTable(name,ascii)) res = false;
  }

  if ( theSubLambdaTable  && isIonisation) {
    const G4String name = 
      GetPhysicsTableFileName(part,directory,"SubLambda",ascii);
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
    G4cout << "Fail to store Physics Tables for " 
	   << particle->GetParticleName()
           << " and process " << GetProcessName()
	   << " in the directory <" << directory
	   << "> " << G4endl;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool G4VEnergyLossProcess::RetrievePhysicsTable(
       const G4ParticleDefinition* part, const G4String& directory,
       G4bool ascii)
{
  G4bool res = true;
  const G4String particleName = part->GetParticleName();

  if(1 < verboseLevel) {
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
      if(yes) yes = G4PhysicsTableHelper::RetrievePhysicsTable(
		    theDEDXTable,filename,ascii);
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "DEDX table for " << particleName 
		 << " is Retrieved from <"
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
      if(yes) yes = G4PhysicsTableHelper::RetrievePhysicsTable(
                    theRangeTableForLoss,filename,ascii);
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "Range table for loss for " << particleName 
		 << " is Retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
        if(fpi) {
          res = false;
	  G4cout << "Range table for loss for " << particleName 
		 << " from file <"
		 << filename << "> is not Retrieved"
		 << G4endl;
        }
      }

      filename = GetPhysicsTableFileName(part,directory,"DEDXnr",ascii);
      yes = theDEDXunRestrictedTable->ExistPhysicsTable(filename);
      if(yes) yes = G4PhysicsTableHelper::RetrievePhysicsTable(
                    theDEDXunRestrictedTable,filename,ascii);
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "Non-restricted DEDX table for " << particleName 
		 << " is Retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
        if (1 < verboseLevel) {
          G4cout << "Non-restricted DEDX table for " << particleName 
		 << " from file <"
                 << filename << "> is not Retrieved"
                 << G4endl;
        }
      }

      filename = GetPhysicsTableFileName(part,directory,"CSDARange",ascii);
      yes = theCSDARangeTable->ExistPhysicsTable(filename);
      if(yes) yes = G4PhysicsTableHelper::RetrievePhysicsTable(
                    theCSDARangeTable,filename,ascii);
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "CSDA Range table for " << particleName 
		 << " is Retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
	G4cout << "CSDA Range table for loss for " << particleName 
	       << " does not exist"
	       << G4endl;
      }

      filename = GetPhysicsTableFileName(part,directory,"InverseRange",ascii);
      yes = theInverseRangeTable->ExistPhysicsTable(filename);
      if(yes)  yes = G4PhysicsTableHelper::RetrievePhysicsTable(
                     theInverseRangeTable,filename,ascii);
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "InverseRange table for " << particleName 
		 << " is Retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
        if(fpi) {
          res = false;
          G4cout << "InverseRange table for " << particleName 
		 << " from file <"
                 << filename << "> is not Retrieved"
                 << G4endl;

        }
      }

      filename = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
      yes = theLambdaTable->ExistPhysicsTable(filename);
      if(yes) yes = G4PhysicsTableHelper::RetrievePhysicsTable(
                    theLambdaTable,filename,ascii);
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "Lambda table for " << particleName 
		 << " is Retrieved from <"
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

      filename = GetPhysicsTableFileName(part,directory,"SubDEDX",ascii);
      yes = theDEDXSubTable->ExistPhysicsTable(filename);
      if(yes) yes = G4PhysicsTableHelper::RetrievePhysicsTable(
                    theDEDXSubTable,filename,ascii);
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "SubDEDX table for " << particleName 
		 << " is Retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
        if(nSCoffRegions) {
          res=false;
          G4cout << "SubDEDX table for " << particleName << " from file <"
                 << filename << "> is not Retrieved"
                 << G4endl;
	}
      }

      filename = GetPhysicsTableFileName(part,directory,"SubLambda",ascii);
      yes = theSubLambdaTable->ExistPhysicsTable(filename);
      if(yes) yes = G4PhysicsTableHelper::RetrievePhysicsTable(
                    theSubLambdaTable,filename,ascii);
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "SubLambda table for " << particleName 
		 << " is Retrieved from <"
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

      filename = GetPhysicsTableFileName(part,directory,"Ionisation",ascii);
      yes = theIonisationTable->ExistPhysicsTable(filename);
      if(yes) {
	theIonisationTable = 
	  G4PhysicsTableHelper::PreparePhysicsTable(theIonisationTable);
        
	yes = G4PhysicsTableHelper::RetrievePhysicsTable(
                    theIonisationTable,filename,ascii);
      }
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "Ionisation table for " << particleName 
		 << " is Retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } 

      filename = GetPhysicsTableFileName(part,directory,"SubIonisation",ascii);
      yes = theIonisationSubTable->ExistPhysicsTable(filename);
      if(yes) {
	theIonisationSubTable = 
	  G4PhysicsTableHelper::PreparePhysicsTable(theIonisationSubTable);
        yes = G4PhysicsTableHelper::RetrievePhysicsTable(
                    theIonisationSubTable,filename,ascii);
      }
      if(yes) {
        if (0 < verboseLevel) {
          G4cout << "SubIonisation table for " << particleName 
		 << " is Retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } 
    }
  }

  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void G4VEnergyLossProcess::AddCollaborativeProcess(
            G4VEnergyLossProcess* p)
{
  G4bool add = true;
  if(nProcesses > 0) {
    for(G4int i=0; i<nProcesses; i++) {
      if(p == scProcesses[i]) {
        add = false;
        break;
      }
    }
  }
  if(add) {
    scProcesses.push_back(p);
    nProcesses++;
    if (0 < verboseLevel) 
      G4cout << "### The process " << p->GetProcessName() 
	     << " is added to the list of collaborative processes of "
	     << GetProcessName() << G4endl; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::GetDEDXDispersion(
                                  const G4MaterialCutsCouple *couple,
                                  const G4DynamicParticle* dp,
                                        G4double length)
{
  DefineMaterial(couple);
  G4double ekin = dp->GetKineticEnergy();
  G4VEmModel* currentModel = SelectModel(ekin*massRatio);
  G4double tmax = currentModel->MaxSecondaryKinEnergy(dp);
  tmax = std::min(tmax,(*theCuts)[currentMaterialIndex]);
  G4double d = 0.0;
  G4VEmFluctuationModel* fm = currentModel->GetModelOfFluctuations();
  if(fm) d = fm->Dispersion(currentMaterial,dp,tmax,length);
  return d;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::ActivateDeexcitation(G4bool, const G4Region*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

