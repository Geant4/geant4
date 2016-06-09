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
// $Id: G4VEnergyLossProcess.cc,v 1.174 2010-12-27 17:42:21 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// 24-06-09 Removed hidden bin in G4PhysicsVector (V.Ivanchenko)
// 15-10-10 Fixed 4-momentum balance if deexcitation is active (L.Pandola)
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
#include "G4EmConfigurator.hh"
#include "G4VAtomDeexcitation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossProcess::G4VEnergyLossProcess(const G4String& name, 
					   G4ProcessType type): 
  G4VContinuousDiscreteProcess(name, type),
  secondaryParticle(0),
  nSCoffRegions(0),
  nDERegions(0),
  idxSCoffRegions(0),
  idxDERegions(0),
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
  baseParticle(0),
  minSubRange(0.1),
  lossFluctuationFlag(true),
  rndmStepFlag(false),
  tablesAreBuilt(false),
  integral(true),
  isIon(false),
  isIonisation(true),
  useSubCutoff(false),
  useDeexcitation(false),
  particle(0),
  currentCouple(0),
  nWarnings(0),
  mfpKinEnergy(0.0)
{
  SetVerboseLevel(1);

  // low energy limit
  lowestKinEnergy  = 1.*eV;

  // Size of tables assuming spline
  minKinEnergy     = 0.1*keV;
  maxKinEnergy     = 10.0*TeV;
  nBins            = 77;
  maxKinEnergyCSDA = 1.0*GeV;
  nBinsCSDA        = 35;

  // default linear loss limit for spline
  linLossLimit  = 0.01;

  // default dRoverRange and finalRange
  SetStepFunction(0.2, 1.0*mm);

  // default lambda factor
  lambdaFactor  = 0.8;

  // particle types
  theElectron   = G4Electron::Electron();
  thePositron   = G4Positron::Positron();
  theGenericIon = 0;

  // run time objects
  pParticleChange = &fParticleChange;
  modelManager = new G4EmModelManager();
  safetyHelper = G4TransportationManager::GetTransportationManager()
    ->GetSafetyHelper();
  aGPILSelection = CandidateForSelection;

  // initialise model
  (G4LossTableManager::Instance())->Register(this);
  fluctModel = 0;
  atomDeexcitation = 0;

  scTracks.reserve(5);
  secParticles.reserve(5);

  // Data for stragling of ranges from ICRU'37 report
  //  const G4int nrbins = 7;
  //vstrag = new G4PhysicsLogVector(keV, GeV, nrbins-1);
  //vstrag->SetSpline(true);
  //G4double s[nrbins] = {-0.2, -0.85, -1.3, -1.578, -1.76, -1.85, -1.9};
  //for(G4int i=0; i<nrbins; ++i) {vstrag->PutValue(i, s[i]);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossProcess::~G4VEnergyLossProcess()
{
  if(1 < verboseLevel) 
    G4cout << "G4VEnergyLossProcess destruct " << GetProcessName() 
	   << G4endl;
  //delete vstrag;
  Clean();

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

void G4VEnergyLossProcess::Clean()
{
  if(1 < verboseLevel) { 
    G4cout << "G4VEnergyLossProcess::Clear() for " << GetProcessName() 
	   << G4endl;
  }
  delete [] theDEDXAtMaxEnergy;
  delete [] theRangeAtMaxEnergy;
  delete [] theEnergyOfCrossSectionMax;
  delete [] theCrossSectionMax;
  delete [] idxSCoffRegions;
  delete [] idxDERegions;

  theDEDXAtMaxEnergy = 0;
  theRangeAtMaxEnergy = 0;
  theEnergyOfCrossSectionMax = 0;
  theCrossSectionMax = 0;
  tablesAreBuilt = false;

  //scTracks.clear();
  scProcesses.clear();
  nProcesses = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::MinPrimaryEnergy(const G4ParticleDefinition*, 
						const G4Material*, 
						G4double cut)
{
  return cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::AddEmModel(G4int order, G4VEmModel* p, 
				      G4VEmFluctuationModel* fluc,
				      const G4Region* region)
{
  modelManager->AddEmModel(order, p, fluc, region);
  if(p) { p->SetParticleChange(pParticleChange, fluc); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::UpdateEmModel(const G4String& nam, 
					 G4double emin, G4double emax)
{
  modelManager->UpdateEmModel(nam, emin, emax);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetEmModel(G4VEmModel* p, G4int index)
{
  G4int n = emModels.size();
  if(index >= n) { for(G4int i=n; i<=index; ++i) {emModels.push_back(0);} }
  emModels[index] = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmModel* G4VEnergyLossProcess::EmModel(G4int index)
{
  G4VEmModel* p = 0;
  if(index >= 0 && index <  G4int(emModels.size())) { p = emModels[index]; }
  return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmModel* G4VEnergyLossProcess::GetModelByIndex(G4int idx, G4bool ver)
{
  return modelManager->GetModel(idx, ver);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4VEnergyLossProcess::NumberOfModels()
{
  return modelManager->NumberOfModels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4VEnergyLossProcess::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::PreparePhysicsTable for "
           << GetProcessName()
           << " for " << part.GetParticleName()
           << G4endl;
  }

  currentCouple = 0;
  preStepLambda = 0.0;
  mfpKinEnergy  = DBL_MAX;
  fRange        = DBL_MAX;
  preStepKinEnergy = 0.0;
  chargeSqRatio = 1.0;
  massRatio = 1.0;
  reduceFactor = 1.0;

  G4LossTableManager* lManager = G4LossTableManager::Instance();

  // Are particle defined?
  if( !particle ) {
    particle = &part;
    if(part.GetParticleType() == "nucleus") {
      if(!theGenericIon) { theGenericIon = G4GenericIon::GenericIon(); }
      if(particle == theGenericIon) { isIon = true; }
      else if(part.GetPDGCharge() > eplus) {
	isIon = true; 

	// generic ions created on-fly
	if(part.GetPDGCharge() > 2.5*eplus) {
	  particle = theGenericIon;
	}
      }
    }
  }

  if( particle != &part) {
    if(part.GetParticleType() == "nucleus") {
      isIon = true;
      lManager->RegisterIon(&part, this);
    } else { 
      lManager->RegisterExtraParticle(&part, this);
    }
    if(1 < verboseLevel) {
      G4cout << "### G4VEnergyLossProcess::PreparePhysicsTable() end for "
	     << part.GetParticleName() << G4endl;
    }
    return;
  }

  Clean();
  lManager->PreparePhysicsTable(&part, this);

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

  if (baseParticle) {
    massRatio = (baseParticle->GetPDGMass())/initialMass;
    G4double q = initialCharge/baseParticle->GetPDGCharge();
    chargeSqRatio = q*q;
    if(chargeSqRatio > 0.0) { reduceFactor = 1.0/(chargeSqRatio*massRatio); }
  }

  // initialisation of models
  G4int nmod = modelManager->NumberOfModels();
  for(G4int i=0; i<nmod; ++i) {
    G4VEmModel* mod = modelManager->GetModel(i);
    if(mod->HighEnergyLimit() > maxKinEnergy) {
      mod->SetHighEnergyLimit(maxKinEnergy);
    }
  }

  theCuts = modelManager->Initialise(particle, secondaryParticle, 
				     minSubRange, verboseLevel);

  // Sub Cutoff and Deexcitation
  if (nSCoffRegions>0 || nDERegions>0) {
    theSubCuts = modelManager->SubCutoff();

    const G4ProductionCutsTable* theCoupleTable=
          G4ProductionCutsTable::GetProductionCutsTable();
    size_t numOfCouples = theCoupleTable->GetTableSize();

    if(nSCoffRegions>0) idxSCoffRegions = new G4bool[numOfCouples];
    if(nDERegions>0) idxDERegions = new G4bool[numOfCouples];
  
    for (size_t j=0; j<numOfCouples; ++j) {

      const G4MaterialCutsCouple* couple = 
	theCoupleTable->GetMaterialCutsCouple(j);
      const G4ProductionCuts* pcuts = couple->GetProductionCuts();
      
      if(nSCoffRegions>0) {
	G4bool reg = false;
	for(G4int i=0; i<nSCoffRegions; ++i) {
	  if( pcuts == scoffRegions[i]->GetProductionCuts()) { reg = true; }
	}
	idxSCoffRegions[j] = reg;
      }
      if(nDERegions>0) {
	G4bool reg = false;
	for(G4int i=0; i<nDERegions; ++i) {
	  if( pcuts == deRegions[i]->GetProductionCuts()) { reg = true; }
	}
	idxDERegions[j] = reg;
      }
    }
  }

  if (1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::Initialise() is done "
           << " for local " << particle->GetParticleName()
	   << " isIon= " << isIon
           << " chargeSqRatio= " << chargeSqRatio
           << " massRatio= " << massRatio
           << " reduceFactor= " << reduceFactor << G4endl;
    if (nSCoffRegions) {
      G4cout << " SubCutoff Regime is ON for regions: " << G4endl;
      for (G4int i=0; i<nSCoffRegions; ++i) {
        const G4Region* r = scoffRegions[i];
	G4cout << "           " << r->GetName() << G4endl;
      }
    }
    if (nDERegions) {
      G4cout << " Deexcitation is ON for regions: " << G4endl;
      for (G4int i=0; i<nDERegions; ++i) {
        const G4Region* r = deRegions[i];
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

  if(&part == particle) {
    if(!tablesAreBuilt) {
      G4LossTableManager::Instance()->BuildPhysicsTable(particle, this);
    }
    if(!baseParticle) {
      if(0 < verboseLevel) PrintInfoDefinition();
    
      // needs to be done only once
      safetyHelper->InitialiseHelper();
    }
  }

  // Added tracking cut to avoid tracking artifacts
  if(isIonisation) { 
    fParticleChange.SetLowEnergyLimit(lowestKinEnergy); 
    atomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();
    if(atomDeexcitation) { 
      if(atomDeexcitation->IsPIXEActive()) { useDeexcitation = true; } 
    }
  }

  if(1 < verboseLevel) {
    G4cout << "### G4VEnergyLossProcess::BuildPhysicsTable() done for "
           << GetProcessName()
           << " and particle " << part.GetParticleName();
    if(isIonisation) { G4cout << "  isIonisation  flag = 1"; }
    G4cout << G4endl;
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
           << " maxKinEnergy= " << emax
           << " nbin= " << bin
           << " EmTableType= " << tType
           << " table= " << table
           << G4endl;
  }
  if(!table) return table;

  G4bool splineFlag = (G4LossTableManager::Instance())->SplineFlag();
  G4PhysicsLogVector* aVector = 0;
  G4PhysicsLogVector* bVector = 0;

  for(size_t i=0; i<numOfCouples; ++i) {

    if(1 < verboseLevel) {
      G4cout << "G4VEnergyLossProcess::BuildDEDXVector flag=  " 
	     << table->GetFlag(i) << G4endl;
    }
    if (table->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = 
	theCoupleTable->GetMaterialCutsCouple(i);
      if(!bVector) {
	aVector = new G4PhysicsLogVector(minKinEnergy, emax, bin);
        bVector = aVector;
      } else {
        aVector = new G4PhysicsLogVector(*bVector);
      }
      aVector->SetSpline(splineFlag);

      modelManager->FillDEDXVector(aVector, couple, tType);
      if(splineFlag) { aVector->FillSecondDerivatives(); }

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
  if(!table) {return table;}

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  G4bool splineFlag = (G4LossTableManager::Instance())->SplineFlag();
  G4PhysicsLogVector* aVector = 0;
  G4double scale = std::log(maxKinEnergy/minKinEnergy);

  for(size_t i=0; i<numOfCouples; ++i) {

    if (table->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = 
	theCoupleTable->GetMaterialCutsCouple(i);
      G4double emin = MinPrimaryEnergy(particle,couple->GetMaterial(),(*theCuts)[i]);
      if(0.0 >= emin) { emin = eV; }
      else if(maxKinEnergy <= emin) { emin = 0.5*maxKinEnergy; }
      G4int bin = G4int(nBins*std::log(maxKinEnergy/emin)/scale + 0.5);
      if(bin < 3) { bin = 3; }
      aVector = new G4PhysicsLogVector(emin, maxKinEnergy, bin);
      aVector->SetSpline(splineFlag);

      modelManager->FillLambdaVector(aVector, couple, true, tType);
      if(splineFlag) { aVector->FillSecondDerivatives(); }

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

void G4VEnergyLossProcess::PrintInfoDefinition()
{
  if(0 < verboseLevel) {
    G4cout << G4endl << GetProcessName() << ":   for  "
           << particle->GetParticleName()
	   << "    SubType= " << GetProcessSubType() 
           << G4endl
           << "      dE/dx and range tables from "
 	   << G4BestUnit(minKinEnergy,"Energy")
           << " to " << G4BestUnit(maxKinEnergy,"Energy")
           << " in " << nBins << " bins" << G4endl
           << "      Lambda tables from threshold to "
           << G4BestUnit(maxKinEnergy,"Energy")
           << " in " << nBins << " bins, spline: " 
	   << (G4LossTableManager::Instance())->SplineFlag()
           << G4endl;
    if(theRangeTableForLoss && isIonisation) {
      G4cout << "      finalRange(mm)= " << finalRange/mm
             << ", dRoverRange= " << dRoverRange
             << ", integral: " << integral
             << ", fluct: " << lossFluctuationFlag
	     << ", linLossLimit= " << linLossLimit
             << G4endl;
    }
    PrintInfo();
    modelManager->DumpModelList(verboseLevel);
    if(theCSDARangeTable && isIonisation) {
      G4cout << "      CSDA range table up"
             << " to " << G4BestUnit(maxKinEnergyCSDA,"Energy")
             << " in " << nBinsCSDA << " bins" << G4endl;
    }
    if(nSCoffRegions>0 && isIonisation) {
      G4cout << "      Subcutoff sampling in " << nSCoffRegions 
	     << " regions" << G4endl;
    }
    if(2 < verboseLevel) {
      G4cout << "      DEDXTable address= " << theDEDXTable << G4endl;
      if(theDEDXTable && isIonisation) G4cout << (*theDEDXTable) << G4endl;
      G4cout << "non restricted DEDXTable address= " 
	     << theDEDXunRestrictedTable << G4endl;
      if(theDEDXunRestrictedTable && isIonisation) {
           G4cout << (*theDEDXunRestrictedTable) << G4endl;
      }
      if(theDEDXSubTable && isIonisation) {
	G4cout << (*theDEDXSubTable) << G4endl;
      }
      G4cout << "      CSDARangeTable address= " << theCSDARangeTable 
	     << G4endl;
      if(theCSDARangeTable && isIonisation) {
	G4cout << (*theCSDARangeTable) << G4endl;
      }
      G4cout << "      RangeTableForLoss address= " << theRangeTableForLoss 
	     << G4endl;
      if(theRangeTableForLoss && isIonisation) {
             G4cout << (*theRangeTableForLoss) << G4endl;
      }
      G4cout << "      InverseRangeTable address= " << theInverseRangeTable 
	     << G4endl;
      if(theInverseRangeTable && isIonisation) {
             G4cout << (*theInverseRangeTable) << G4endl;
      }
      G4cout << "      LambdaTable address= " << theLambdaTable << G4endl;
      if(theLambdaTable && isIonisation) {
	G4cout << (*theLambdaTable) << G4endl;
      }
      G4cout << "      SubLambdaTable address= " << theSubLambdaTable << G4endl;
      if(theSubLambdaTable && isIonisation) {
	G4cout << (*theSubLambdaTable) << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::ActivateSubCutoff(G4bool val, const G4Region* r)
{
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  const G4Region* reg = r;
  if (!reg) {reg = regionStore->GetRegion("DefaultRegionForTheWorld", false);}

  // the region is in the list
  if (nSCoffRegions) {
    for (G4int i=0; i<nSCoffRegions; ++i) {
      if (reg == scoffRegions[i]) {
	if(!val) deRegions[i] = 0;
        return;
      }
    }
  }

  // new region 
  if(val) {
    useSubCutoff = true;
    scoffRegions.push_back(reg);
    ++nSCoffRegions;
  } else {
    useSubCutoff = false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::ActivateDeexcitation(G4bool val, const G4Region* r)
{
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  const G4Region* reg = r;
  if (!reg) {reg = regionStore->GetRegion("DefaultRegionForTheWorld", false);}

  // the region is in the list
  if (nDERegions) {
    for (G4int i=0; i<nDERegions; ++i) {
      if (reg == deRegions[i]) {
	if(!val) { deRegions[i] = 0; }
        return;
      }
    }
  }

  // new region 
  if(val) {
    useDeexcitation = true;
    deRegions.push_back(reg);
    ++nDERegions;
  } else {
    useDeexcitation = false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::AlongStepGetPhysicalInteractionLength(
                             const G4Track&,G4double,G4double,G4double&,
                             G4GPILSelection* selection)
{
  G4double x = DBL_MAX;
  *selection = aGPILSelection;
  if(isIonisation) {
    fRange = GetScaledRangeForScaledEnergy(preStepScaledEnergy)*reduceFactor;

    x = fRange;
    G4double y = x*dRoverRange;
    G4double finR = finalRange;
    if(rndmStepFlag) { 
      finR = std::min(finR,currentCouple->GetProductionCuts()->GetProductionCut(1)); 
    }
    if(x > finR) { x = y + finR*(1.0 - dRoverRange)*(2.0 - finR/fRange); }
    //if(x > finalRange && y < currentMinStep) { 
    //  x = y + finalRange*(1.0 - dRoverRange)*(2.0 - finalRange/fRange);
    //} else if (rndmStepFlag) { x = SampleRange(); }
    //G4cout<<GetProcessName()<<": e= "<<preStepKinEnergy
    //   <<" range= "<<fRange <<" cMinSt="<<currentMinStep
    //	  << " y= " << y << " finR= " << prec
    //    << " limit= " << x <<G4endl;
  }
  //G4cout<<GetProcessName()<<": e= "<<preStepKinEnergy
  //  <<" stepLimit= "<<x<<G4endl;
  return x;
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

  // initialisation of material, mass, charge, model at the beginning of the step
  DefineMaterial(track.GetMaterialCutsCouple());

  const G4ParticleDefinition* currPart = track.GetParticleDefinition();
  if(theGenericIon == particle) {
    massRatio = proton_mass_c2/currPart->GetPDGMass();
  }  
  preStepKinEnergy    = track.GetKineticEnergy();
  preStepScaledEnergy = preStepKinEnergy*massRatio;
  SelectModel(preStepScaledEnergy);
  if(!currentModel->IsActive(preStepScaledEnergy)) { return x; }

  if(isIon) { 
    chargeSqRatio = currentModel->ChargeSquareRatio(track);
    reduceFactor  = 1.0/(chargeSqRatio*massRatio);
  }
  //G4cout << "q2= " << chargeSqRatio << " massRatio= " << massRatio << G4endl; 
  // initialisation for sampling of the interaction length 
  if(previousStepSize <= DBL_MIN) { theNumberOfInteractionLengthLeft = -1.0; }
  if(theNumberOfInteractionLengthLeft < 0.0) { mfpKinEnergy = DBL_MAX; }

  // compute mean free path
  if(preStepScaledEnergy < mfpKinEnergy) {
    if (integral) { ComputeLambdaForScaledEnergy(preStepScaledEnergy); }
    else  { preStepLambda = GetLambdaForScaledEnergy(preStepScaledEnergy); }
    if(preStepLambda <= DBL_MIN) { mfpKinEnergy = 0.0; }
  }

  // non-zero cross section
  if(preStepLambda > DBL_MIN) { 
    if (theNumberOfInteractionLengthLeft < 0.0) {
      // beggining of tracking (or just after DoIt of this process)
      //G4cout<<"G4VEnergyLossProcess::PostStepGetPhysicalInteractionLength Reset"<<G4endl;
      ResetNumberOfInteractionLengthLeft();
    } else if(currentInteractionLength < DBL_MAX) {
      // subtract NumberOfInteractionLengthLeft
      SubtractNumberOfInteractionLengthLeft(previousStepSize);
      if(theNumberOfInteractionLengthLeft < 0.) {
	theNumberOfInteractionLengthLeft = perMillion;
      }
    }

    // get mean free path and step limit
    currentInteractionLength = 1.0/preStepLambda;
    x = theNumberOfInteractionLengthLeft * currentInteractionLength;

#ifdef G4VERBOSE
    if (verboseLevel>2){
      G4cout << "G4VEnergyLossProcess::PostStepGetPhysicalInteractionLength ";
      G4cout << "[ " << GetProcessName() << "]" << G4endl; 
      G4cout << " for " << currPart->GetParticleName() 
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
      if(theNumberOfInteractionLengthLeft < 0.) {
	theNumberOfInteractionLengthLeft = perMillion;
      }
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
  if(!isIonisation || !currentModel->IsActive(preStepScaledEnergy)) {
    return &fParticleChange;
  }

  // Get the actual (true) Step length
  G4double length = step.GetStepLength();
  if(length <= DBL_MIN) { return &fParticleChange; }
  G4double eloss  = 0.0;
 
  /*  
  if(-1 < verboseLevel) {
    const G4ParticleDefinition* d = track.GetParticleDefinition();
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

  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();

  // stopping
  if (length >= fRange) {
    eloss = preStepKinEnergy;
    if (useDeexcitation) {
      if(atomDeexcitation) { 
	atomDeexcitation->AlongStepDeexcitation(&fParticleChange, step, 
						eloss, currentMaterialIndex);
      } else if(idxDERegions[currentMaterialIndex]) {
	currentModel->SampleDeexcitationAlongStep(currentMaterial, track, eloss);
      }
      if(eloss < 0.0) { eloss = 0.0; }
    }
    fParticleChange.SetProposedKineticEnergy(0.0);
    fParticleChange.ProposeLocalEnergyDeposit(eloss);
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
	     << " lim(MeV)= " << preStepKinEnergy*linLossLimit/MeV
             << G4endl;
    */
  }

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

  // SubCutOff 
  if(useSubCutoff) {
    if(idxSCoffRegions[currentMaterialIndex]) {

      G4bool yes = false;
      G4StepPoint* prePoint  = step.GetPreStepPoint();

      // Check boundary
      if(prePoint->GetStepStatus() == fGeomBoundary) { yes = true; }

      // Check PrePoint
      else {
	G4double preSafety  = prePoint->GetSafety();
	G4double rcut = currentCouple->GetProductionCuts()->GetProductionCut(1);

	// recompute presafety
        if(preSafety < rcut) {
	  preSafety = safetyHelper->ComputeSafety(prePoint->GetPosition());
	}

        if(preSafety < rcut) { yes = true; }

	// Check PostPoint
	else {
	  G4double postSafety = preSafety - length; 
	  if(postSafety < rcut) {
	    postSafety = 
	      safetyHelper->ComputeSafety(step.GetPostStepPoint()->GetPosition());
	    if(postSafety < rcut) { yes = true; }
	  }
	}
      }
  
      // Decide to start subcut sampling
      if(yes) {

        cut = (*theSubCuts)[currentMaterialIndex];
 	eloss -= GetSubDEDXForScaledEnergy(preStepScaledEnergy)*length;
	scTracks.clear();
	SampleSubCutSecondaries(scTracks, step, 
				currentModel,currentMaterialIndex);
	// add bremsstrahlung sampling
	/*
	if(nProcesses > 0) {
	  for(G4int i=0; i<nProcesses; ++i) {
	    (scProcesses[i])->SampleSubCutSecondaries(
		scTracks, step, (scProcesses[i])->
		SelectModelForMaterial(preStepKinEnergy, currentMaterialIndex),
		currentMaterialIndex);
	  }
	} 
	*/   
	G4int n = scTracks.size();
	if(n>0) {
	  G4ThreeVector mom = dynParticle->GetMomentum();
	  fParticleChange.SetNumberOfSecondaries(n);
	  for(G4int i=0; i<n; ++i) {
	    G4Track* t = scTracks[i];
	    G4double e = t->GetKineticEnergy();
	    if (t->GetParticleDefinition() == thePositron) { 
	      e += 2.0*electron_mass_c2; 
	    }
	    esec += e;
	    pParticleChange->AddSecondary(t);
	  }      
	}
      }
    }
  }

  // Corrections, which cannot be tabulated
  if(isIon) {
    G4double eadd = 0.0;
    currentModel->CorrectionsAlongStep(currentCouple, dynParticle, 
				       eloss, eadd, length);
  }

  // Sample fluctuations
  if (lossFluctuationFlag) {
    G4VEmFluctuationModel* fluc = currentModel->GetModelOfFluctuations();
    if(fluc && 
      (eloss + esec + lowestKinEnergy) < preStepKinEnergy) {

      G4double tmax = 
	std::min(currentModel->MaxSecondaryKinEnergy(dynParticle),cut);
      G4double emean = eloss;
      eloss = fluc->SampleFluctuations(currentMaterial,dynParticle,
				       tmax,length,emean);
      /*                            
      if(-1 < verboseLevel) 
      G4cout << "After fluct: eloss(MeV)= " << eloss/MeV
             << " fluc= " << (eloss-eloss0)/MeV
             << " ChargeSqRatio= " << chargeSqRatio
             << " massRatio= " << massRatio
             << " tmax= " << tmax
             << G4endl;
      */
    }
  }

  // deexcitation
  if (useDeexcitation) {
    G4double eloss_before = eloss;
    if(atomDeexcitation) { 
      atomDeexcitation->AlongStepDeexcitation(&fParticleChange, step, 
					      eloss, currentMaterialIndex);
    } else if(idxDERegions[currentMaterialIndex] && eloss > 0.0) {
      currentModel->SampleDeexcitationAlongStep(currentMaterial, track, eloss);
    }
    esec += eloss_before - eloss;
  }

  // Energy balanse
  G4double finalT = preStepKinEnergy - eloss - esec;
  if (finalT <= lowestKinEnergy) {
    eloss += finalT;
    finalT = 0.0;
  } else if(isIon) {
    fParticleChange.SetProposedCharge(
      currentModel->GetParticleCharge(track.GetParticleDefinition(),
				      currentMaterial,finalT));
  }

  if(eloss < 0.0) { eloss = 0.0; }
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
       G4int idx) 
{
  // Fast check weather subcutoff can work
  G4double subcut = (*theSubCuts)[idx];
  G4double cut = (*theCuts)[idx];
  if(cut <= subcut) { return; }

  const G4Track* track = step.GetTrack();
  const G4DynamicParticle* dp = track->GetDynamicParticle();
  G4double e = dp->GetKineticEnergy()*massRatio;
  G4double cross = chargeSqRatio*(((*theSubLambdaTable)[idx])->Value(e));
  G4double length = step.GetStepLength();

  // negligible probability to get any interaction
  if(length*cross < perMillion) { return; }
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
  G4StepPoint* postStepPoint = step.GetPostStepPoint();
  G4ThreeVector prepoint = preStepPoint->GetPosition();
  G4ThreeVector dr = postStepPoint->GetPosition() - prepoint;
  G4double pretime = preStepPoint->GetGlobalTime();
  G4double dt = postStepPoint->GetGlobalTime() - pretime;
  //G4double dt = length/preStepPoint->GetVelocity();
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
    for(it=secParticles.begin(); it!=secParticles.end(); ++it) {

      G4bool addSec = true;
      /*
      // do not track very low-energy delta-electrons
      if(theSecondaryRangeTable && (*it)->GetParticleDefinition() == theElectron) {
	G4double ekin = (*it)->GetKineticEnergy();
	G4double rg = ((*theSecondaryRangeTable)[idx]->Value(ekin));
	//          if(rg < currentMinSafety) {
	if(rg < safetyHelper->ComputeSafety(r)) {
	  extraEdep += ekin;
	  delete (*it);
	  addSec = false;
	}
      }
      */
      if(addSec) {
	G4Track* t = new G4Track((*it), pretime + fragment*dt, r);
	t->SetTouchableHandle(track->GetTouchableHandle());
	tracks.push_back(t);

	/*	
	if(-1 < verboseLevel) 
	  G4cout << "New track " << t->GetParticleDefinition()->GetParticleName()
		 << " e(keV)= " << t->GetKineticEnergy()/keV
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
  // In all cases clear number of interaction lengths
  theNumberOfInteractionLengthLeft = -1.0;

  fParticleChange.InitializeForPostStep(track);
  G4double finalT = track.GetKineticEnergy();
  if(finalT <= lowestKinEnergy) { return &fParticleChange; }

  G4double postStepScaledEnergy = finalT*massRatio;

  if(!currentModel->IsActive(postStepScaledEnergy)) { return &fParticleChange; }
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
      ++nWarnings;
    }
    */
    if(lx <= 0.0) {
      return &fParticleChange;
    } else if(preStepLambda*G4UniformRand() > lx) {
      return &fParticleChange;
    }
  }

  SelectModel(postStepScaledEnergy);
  if(useDeexcitation && !atomDeexcitation) { 
    currentModel->SetDeexcitationFlag(idxDERegions[currentMaterialIndex]);
  }

  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  G4double tcut = (*theCuts)[currentMaterialIndex];

  // sample secondaries
  secParticles.clear();
  currentModel->SampleSecondaries(&secParticles, currentCouple, dynParticle, tcut);

  // save secondaries
  G4int num = secParticles.size();
  if(num > 0) {
    fParticleChange.SetNumberOfSecondaries(num);
    for (G4int i=0; i<num; ++i) {
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
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4VEnergyLossProcess::StorePhysicsTable(
       const G4ParticleDefinition* part, const G4String& directory, 
       G4bool ascii)
{
  G4bool res = true;
  if ( baseParticle || part != particle ) return res;

  if(!StoreTable(part,theDEDXTable,ascii,directory,"DEDX")) 
    {res = false;}

  if(!StoreTable(part,theDEDXunRestrictedTable,ascii,directory,"DEDXnr")) 
    {res = false;}

  if(!StoreTable(part,theDEDXSubTable,ascii,directory,"SubDEDX")) 
    {res = false;}

  if(!StoreTable(part,theIonisationTable,ascii,directory,"Ionisation")) 
    {res = false;}

  if(!StoreTable(part,theIonisationSubTable,ascii,directory,"SubIonisation")) 
    {res = false;}

  if(isIonisation &&
     !StoreTable(part,theCSDARangeTable,ascii,directory,"CSDARange")) 
    {res = false;}

  if(isIonisation &&
     !StoreTable(part,theRangeTableForLoss,ascii,directory,"Range")) 
    {res = false;}
  
  if(isIonisation &&
     !StoreTable(part,theInverseRangeTable,ascii,directory,"InverseRange")) 
    {res = false;}
  
  if(!StoreTable(part,theLambdaTable,ascii,directory,"Lambda")) 
    {res = false;}

  if(!StoreTable(part,theSubLambdaTable,ascii,directory,"SubLambda")) 
    {res = false;}

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

G4bool 
G4VEnergyLossProcess::RetrievePhysicsTable(const G4ParticleDefinition* part, 
					   const G4String& directory,
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

    if ( !baseParticle ) {

      G4bool fpi = true;
      if(!RetrieveTable(part,theDEDXTable,ascii,directory,"DEDX",fpi)) 
	{fpi = false;}

      // ionisation table keeps individual dEdx and not sum of sub-processes
      if(!RetrieveTable(part,theDEDXTable,ascii,directory,"Ionisation",false)) 
	{fpi = false;}

      if(!RetrieveTable(part,theRangeTableForLoss,ascii,directory,"Range",fpi)) 
        {res = false;}

      if(!RetrieveTable(part,theDEDXunRestrictedTable,ascii,directory,"DEDXnr",false)) 
	{res = false;}

      if(!RetrieveTable(part,theCSDARangeTable,ascii,directory,"CSDARange",false)) 
	{res = false;}

      if(!RetrieveTable(part,theInverseRangeTable,ascii,directory,"InverseRange",fpi)) 
        {res = false;}

      if(!RetrieveTable(part,theLambdaTable,ascii,directory,"Lambda",true)) 
        {res = false;}

      G4bool yes = false;
      if(nSCoffRegions > 0) {yes = true;}

      if(!RetrieveTable(part,theDEDXSubTable,ascii,directory,"SubDEDX",yes)) 
        {res = false;}

      if(!RetrieveTable(part,theSubLambdaTable,ascii,directory,"SubLambda",yes)) 
        {res = false;}

      if(!fpi) yes = false;
      if(!RetrieveTable(part,theIonisationSubTable,ascii,directory,"SubIonisation",yes))
        {res = false;}
    }
  }

  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool G4VEnergyLossProcess::StoreTable(const G4ParticleDefinition* part, 
					G4PhysicsTable* aTable, G4bool ascii,
					const G4String& directory,
					const G4String& tname)
{
  G4bool res = true;
  if ( aTable ) {
    const G4String name = GetPhysicsTableFileName(part,directory,tname,ascii);
    if( !aTable->StorePhysicsTable(name,ascii)) res = false;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool 
G4VEnergyLossProcess::RetrieveTable(const G4ParticleDefinition* part, 
				    G4PhysicsTable* aTable, 
				    G4bool ascii,
				    const G4String& directory,
				    const G4String& tname,
				    G4bool mandatory)
{
  G4bool isRetrieved = false;
  G4String filename = GetPhysicsTableFileName(part,directory,tname,ascii);
  if(aTable) {
    if(aTable->ExistPhysicsTable(filename)) {
      if(G4PhysicsTableHelper::RetrievePhysicsTable(aTable,filename,ascii)) {
	isRetrieved = true;
	if((G4LossTableManager::Instance())->SplineFlag()) {
	  size_t n = aTable->length();
	  for(size_t i=0; i<n; ++i) {
	    if((*aTable)[i]) { (*aTable)[i]->SetSpline(true); }
	  }
	}
	if (0 < verboseLevel) {
	  G4cout << tname << " table for " << part->GetParticleName() 
		 << " is Retrieved from <" << filename << ">"
		 << G4endl;
	}
      }
    }
  }
  if(mandatory && !isRetrieved) {
    if(0 < verboseLevel) {
      G4cout << tname << " table for " << part->GetParticleName() 
	     << " from file <"
	     << filename << "> is not Retrieved"
	     << G4endl;
    }
    return false;
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::GetDEDXDispersion(
                                  const G4MaterialCutsCouple *couple,
                                  const G4DynamicParticle* dp,
                                        G4double length)
{
  DefineMaterial(couple);
  G4double ekin = dp->GetKineticEnergy();
  SelectModel(ekin*massRatio);
  G4double tmax = currentModel->MaxSecondaryKinEnergy(dp);
  tmax = std::min(tmax,(*theCuts)[currentMaterialIndex]);
  G4double d = 0.0;
  G4VEmFluctuationModel* fm = currentModel->GetModelOfFluctuations();
  if(fm) d = fm->Dispersion(currentMaterial,dp,tmax,length);
  return d;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::CrossSectionPerVolume(
	 G4double kineticEnergy, const G4MaterialCutsCouple* couple)
{
  // Cross section per volume is calculated
  DefineMaterial(couple);
  G4double cross = 0.0;
  if(theLambdaTable) {
    cross = ((*theLambdaTable)[currentMaterialIndex])->Value(kineticEnergy);
  } else {
    SelectModel(kineticEnergy);
    cross = currentModel->CrossSectionPerVolume(currentMaterial,
						particle, kineticEnergy,
						(*theCuts)[currentMaterialIndex]);
  }
  if(cross < 0.0) { cross = 0.0; }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::MeanFreePath(const G4Track& track)
{
  DefineMaterial(track.GetMaterialCutsCouple());
  preStepLambda = GetLambdaForScaledEnergy(track.GetKineticEnergy()*massRatio);
  G4double x = DBL_MAX;
  if(DBL_MIN < preStepLambda) { x = 1.0/preStepLambda; }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::ContinuousStepLimit(const G4Track& track, 
						   G4double x, G4double y, 
						   G4double& z)
{
  G4GPILSelection sel;
  return AlongStepGetPhysicalInteractionLength(track, x, y, z, &sel);
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

G4double G4VEnergyLossProcess::GetContinuousStepLimit(
		const G4Track&,
                G4double, G4double, G4double&)
{
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* 
G4VEnergyLossProcess::LambdaPhysicsVector(const G4MaterialCutsCouple* /*couple*/, 
					  G4double /*cut*/)
{
  G4PhysicsVector* v = new G4PhysicsLogVector(minKinEnergy, maxKinEnergy, nBins);
  v->SetSpline((G4LossTableManager::Instance())->SplineFlag());
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void G4VEnergyLossProcess::AddCollaborativeProcess(
            G4VEnergyLossProcess* p)
{
  G4bool add = true;
  if(p->GetProcessName() != "eBrem") { add = false; }
  if(add && nProcesses > 0) {
    for(G4int i=0; i<nProcesses; ++i) {
      if(p == scProcesses[i]) {
        add = false;
        break;
      }
    }
  }
  if(add) {
    scProcesses.push_back(p);
    ++nProcesses;
    if (1 < verboseLevel) { 
      G4cout << "### The process " << p->GetProcessName() 
	     << " is added to the list of collaborative processes of "
	     << GetProcessName() << G4endl; 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetDEDXTable(G4PhysicsTable* p, G4EmTableType tType)
{
  if(fTotal == tType && theDEDXunRestrictedTable != p && !baseParticle) {
    if(theDEDXunRestrictedTable) {
      theDEDXunRestrictedTable->clearAndDestroy();
      delete theDEDXunRestrictedTable;
    } 
    theDEDXunRestrictedTable = p;
    if(p) {
      size_t n = p->length();
      G4PhysicsVector* pv = (*p)[0];
      G4double emax = maxKinEnergyCSDA;
      theDEDXAtMaxEnergy = new G4double [n];

      for (size_t i=0; i<n; ++i) {
	pv = (*p)[i];
	G4double dedx = pv->Value(emax);
	theDEDXAtMaxEnergy[i] = dedx;
	//G4cout << "i= " << i << " emax(MeV)= " << emax/MeV<< " dedx= " 
	//<< dedx << G4endl;
      }
    }

  } else if(fRestricted == tType && theDEDXTable != p) {
    //G4cout << "G4VEnergyLossProcess::SetDEDXTable " << particle->GetParticleName()
    //	   << " old table " << theDEDXTable << " new table " << p 
    //	   << " ion " << theIonisationTable << " bp " << baseParticle << G4endl;
    if(theDEDXTable && !baseParticle) {
      if(theDEDXTable == theIonisationTable) { theIonisationTable = 0; }
      theDEDXTable->clearAndDestroy();
      delete theDEDXTable;
    }
    theDEDXTable = p;
  } else if(fSubRestricted == tType && theDEDXSubTable != p) {    
    if(theDEDXSubTable && !baseParticle) {
      if(theDEDXSubTable == theIonisationSubTable) { theIonisationSubTable = 0; }
      theDEDXSubTable->clearAndDestroy();
      delete theDEDXSubTable;
    }
    theDEDXSubTable = p;
  } else if(fIsIonisation == tType && theIonisationTable != p) {    
    if(theIonisationTable && theIonisationTable != theDEDXTable && !baseParticle) {
      theIonisationTable->clearAndDestroy();
      delete theIonisationTable;
    }
    theIonisationTable = p;
  } else if(fIsSubIonisation == tType && theIonisationSubTable != p) {    
    if(theIonisationSubTable && theIonisationSubTable != theDEDXSubTable && !baseParticle) {
      theIonisationSubTable->clearAndDestroy();
      delete theIonisationSubTable;
    }
    theIonisationSubTable = p;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetCSDARangeTable(G4PhysicsTable* p)
{
  if(theCSDARangeTable != p) { theCSDARangeTable = p; }

  if(p) {
    size_t n = p->length();
    G4PhysicsVector* pv = (*p)[0];
    G4double emax = maxKinEnergyCSDA;
    theRangeAtMaxEnergy = new G4double [n];

    for (size_t i=0; i<n; ++i) {
      pv = (*p)[i];
      G4double r2 = pv->Value(emax);
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
  if(theLambdaTable != p) { theLambdaTable = p; }
  tablesAreBuilt = true;

  if(p) {
    size_t n = p->length();
    G4PhysicsVector* pv = (*p)[0];
    G4double e, s, smax, emax;
    theEnergyOfCrossSectionMax = new G4double [n];
    theCrossSectionMax = new G4double [n];

    for (size_t i=0; i<n; ++i) {
      pv = (*p)[i];
      emax = DBL_MAX;
      smax = 0.0;
      if(pv) {
        size_t nb = pv->GetVectorLength();
        if(nb > 0) {
	  for (size_t j=0; j<nb; ++j) {
	    e = pv->Energy(j);
	    s = (*pv)(j);
	    if(s > smax) {
	      smax = s;
	      emax = e;
	    }
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

const G4Element* G4VEnergyLossProcess::GetCurrentElement() const
{
  const G4Element* elm = 0;
  if(currentModel) { elm = currentModel->GetCurrentElement(); }
  return elm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

