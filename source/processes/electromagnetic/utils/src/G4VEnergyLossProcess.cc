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
// $Id: G4VEnergyLossProcess.cc 104349 2017-05-26 07:18:59Z gcosmo $
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
// 06-03-03 Control on GenericIons using SubType+ update verbose (V.Ivanchenko)
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
// 01-25-09 (Xin Dong) Phase II change for Geant4 multi-threading:
//          New methods SlavePreparePhysicsTable, SlaveBuildPhysicsTable
//          Worker threads share physics tables with the master thread for
//          this kind of process. This member function is used by worker
//          threads to achieve the partial effect of the master thread when
//          it builds physcis tables.
// 15-10-10 Fixed 4-momentum balance if deexcitation is active (L.Pandola)
// 30-05-12 Call new ApplySecondaryBiasing so 2ries may be unique (D. Sawkey)
// 30-05-12 Fix bug in forced biasing: now called on first step (D. Sawkey)
// 04-06-13 Adoptation to MT mode, adding internal cache to GetRangeForLoss,
//          more accurate initialisation for ions (V.Ivanchenko)
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
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4LossTableManager.hh"
#include "G4LossTableBuilder.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4DataVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4VParticleChange.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"
#include "G4UnitsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4SafetyHelper.hh"
#include "G4TransportationManager.hh"
#include "G4EmConfigurator.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4VSubCutProducer.hh"
#include "G4EmBiasingManager.hh"
#include "G4Log.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossProcess::G4VEnergyLossProcess(const G4String& name, 
                                           G4ProcessType type): 
  G4VContinuousDiscreteProcess(name, type),
  secondaryParticle(nullptr),
  nSCoffRegions(0),
  idxSCoffRegions(nullptr),
  nProcesses(0),
  theDEDXTable(nullptr),
  theDEDXSubTable(nullptr),
  theDEDXunRestrictedTable(nullptr),
  theIonisationTable(nullptr),
  theIonisationSubTable(nullptr),
  theRangeTableForLoss(nullptr),
  theCSDARangeTable(nullptr),
  theSecondaryRangeTable(nullptr),
  theInverseRangeTable(nullptr),
  theLambdaTable(nullptr),
  theSubLambdaTable(nullptr),
  theDensityFactor(nullptr),
  theDensityIdx(nullptr),
  baseParticle(nullptr),
  lossFluctuationFlag(true),
  rndmStepFlag(false),
  tablesAreBuilt(false),
  integral(true),
  isIon(false),
  isIonisation(true),
  useSubCutoff(false),
  useDeexcitation(false),
  particle(nullptr),
  currentCouple(nullptr),
  mfpKinEnergy(0.0)
{
  theParameters = G4EmParameters::Instance();
  SetVerboseLevel(1);

  // low energy limit
  lowestKinEnergy  = theParameters->LowestElectronEnergy();
  preStepKinEnergy = 0.0;
  preStepRangeEnergy = 0.0;
  computedRange = DBL_MAX;

  // Size of tables assuming spline
  minKinEnergy     = 0.1*keV;
  maxKinEnergy     = 100.0*TeV;
  nBins            = 77;
  maxKinEnergyCSDA = 1.0*GeV;
  nBinsCSDA        = 35;
  actMinKinEnergy = actMaxKinEnergy = actBinning = actLinLossLimit 
    = actLossFluc = actIntegral = actStepFunc = false;

  // default linear loss limit for spline
  linLossLimit  = 0.01;
  dRoverRange = 0.2;
  finalRange = CLHEP::mm;

  // default lambda factor
  lambdaFactor  = 0.8;

  // cross section biasing
  biasFactor = 1.0;

  // particle types
  theElectron   = G4Electron::Electron();
  thePositron   = G4Positron::Positron();
  theGamma      = G4Gamma::Gamma();
  theGenericIon = nullptr;

  // run time objects
  pParticleChange = &fParticleChange;
  fParticleChange.SetSecondaryWeightByProcess(true);
  modelManager = new G4EmModelManager();
  safetyHelper = G4TransportationManager::GetTransportationManager()
    ->GetSafetyHelper();
  aGPILSelection = CandidateForSelection;

  // initialise model
  lManager = G4LossTableManager::Instance();
  lManager->Register(this);
  fluctModel = nullptr;
  currentModel = nullptr;
  atomDeexcitation = nullptr;
  subcutProducer = nullptr;

  biasManager  = nullptr;
  biasFlag     = false; 
  weightFlag   = false; 
  isMaster     = true;
  lastIdx      = 0;

  idxDEDX = idxDEDXSub = idxDEDXunRestricted = idxIonisation =
    idxIonisationSub = idxRange = idxCSDA = idxSecRange =
    idxInverseRange = idxLambda = idxSubLambda = 0;

  scTracks.reserve(5);
  secParticles.reserve(5);

  theCuts = theSubCuts = nullptr;
  currentMaterial = nullptr;
  currentCoupleIndex  = basedCoupleIndex = 0;
  massRatio = fFactor = reduceFactor = chargeSqRatio = 1.0;
  preStepLambda = preStepScaledEnergy = fRange = 0.0;

  secID = biasID = subsecID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossProcess::~G4VEnergyLossProcess()
{
  /*
  G4cout << "** G4VEnergyLossProcess::~G4VEnergyLossProcess() for " 
         << GetProcessName() << " isMaster: " << isMaster
         << "  basePart: " << baseParticle 
         << G4endl;
  */
  Clean();

  // G4cout << " isIonisation " << isIonisation << "  " 
  //   << theDEDXTable << "  " <<  theIonisationTable << G4endl;

  if (isMaster && !baseParticle) {
    if(theDEDXTable) {

      //G4cout << " theIonisationTable " << theIonisationTable << G4endl;
      if(theIonisationTable == theDEDXTable) { theIonisationTable = 0; }
      //G4cout << " delete theDEDXTable " << theDEDXTable << G4endl;
      theDEDXTable->clearAndDestroy();
      delete theDEDXTable;
      theDEDXTable = nullptr;
      if(theDEDXSubTable) {
        if(theIonisationSubTable == theDEDXSubTable) 
          { theIonisationSubTable = nullptr; }
        theDEDXSubTable->clearAndDestroy();
        delete theDEDXSubTable;
        theDEDXSubTable = nullptr;
      }
    }
    //G4cout << " theIonisationTable " << theIonisationTable << G4endl;
    if(theIonisationTable) {
      //G4cout << " delete theIonisationTable " << theIonisationTable << G4endl;
      theIonisationTable->clearAndDestroy();
      delete theIonisationTable;
      theIonisationTable = nullptr;
    }
    if(theIonisationSubTable) {
      theIonisationSubTable->clearAndDestroy();
      delete theIonisationSubTable;
      theIonisationSubTable = nullptr;
    }
    if(theDEDXunRestrictedTable && isIonisation) {
      theDEDXunRestrictedTable->clearAndDestroy();
      delete theDEDXunRestrictedTable;
      theDEDXunRestrictedTable = nullptr;
    }
    if(theCSDARangeTable && isIonisation) {
      theCSDARangeTable->clearAndDestroy();
      delete theCSDARangeTable;
      theCSDARangeTable = nullptr;
    }
    //G4cout << "delete RangeTable: " << theRangeTableForLoss << G4endl;
    if(theRangeTableForLoss && isIonisation) {
      theRangeTableForLoss->clearAndDestroy();
      delete theRangeTableForLoss;
      theRangeTableForLoss = nullptr;
    }
    //G4cout << "delete InvRangeTable: " << theInverseRangeTable << G4endl;
    if(theInverseRangeTable && isIonisation /*&& !isIon*/) {
      theInverseRangeTable->clearAndDestroy();
      delete theInverseRangeTable;
      theInverseRangeTable = nullptr;
    }
    //G4cout << "delete LambdaTable: " << theLambdaTable << G4endl;
    if(theLambdaTable) {
      theLambdaTable->clearAndDestroy();
      delete theLambdaTable;
      theLambdaTable = nullptr;
    }
    if(theSubLambdaTable) {
      theSubLambdaTable->clearAndDestroy();
      delete theSubLambdaTable;
      theSubLambdaTable = nullptr;
    }
  }
 
  delete modelManager;
  delete biasManager;
  lManager->DeRegister(this);
  //G4cout << "** all removed" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::Clean()
{
  /*
  if(1 < verboseLevel) { 
    G4cout << "G4VEnergyLossProcess::Clear() for " << GetProcessName() 
           << G4endl;
  }
  */
  delete [] idxSCoffRegions;

  tablesAreBuilt = false;

  scProcesses.clear();
  nProcesses = 0;

  idxDEDX = idxDEDXSub = idxDEDXunRestricted = idxIonisation =
    idxIonisationSub = idxRange = idxCSDA = idxSecRange =
    idxInverseRange = idxLambda = idxSubLambda = 0;
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

G4VEmModel* G4VEnergyLossProcess::EmModel(G4int index) const
{
  G4VEmModel* p = nullptr;
  if(index >= 0 && index <  G4int(emModels.size())) { p = emModels[index]; }
  return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmModel* G4VEnergyLossProcess::GetModelByIndex(G4int idx, G4bool ver) const
{
  return modelManager->GetModel(idx, ver);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4VEnergyLossProcess::NumberOfModels() const
{
  return modelManager->NumberOfModels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4VEnergyLossProcess::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::PreparePhysicsTable for "
           << GetProcessName() << " for " << part.GetParticleName() 
           << "  " << this << G4endl;
  }

  const G4VEnergyLossProcess* masterProcess = 
    static_cast<const G4VEnergyLossProcess*>(GetMasterProcess());
  if(masterProcess && masterProcess != this) { isMaster = false; }

  currentCouple = nullptr;
  preStepLambda = 0.0;
  mfpKinEnergy  = DBL_MAX;
  fRange        = DBL_MAX;
  preStepKinEnergy = 0.0;
  preStepRangeEnergy = 0.0;
  chargeSqRatio = 1.0;
  massRatio = 1.0;
  reduceFactor = 1.0;
  fFactor = 1.0;
  lastIdx = 0;

  // Are particle defined?
  if( !particle ) { particle = &part; }

  if(part.GetParticleType() == "nucleus") {

    G4String pname = part.GetParticleName();
    if(pname != "deuteron" && pname != "triton" &&
       pname != "alpha+"   && pname != "helium" &&
       pname != "hydrogen") {

      if(!theGenericIon) {
        theGenericIon = 
          G4ParticleTable::GetParticleTable()->FindParticle("GenericIon");
      }
      isIon = true; 
      if(theGenericIon && particle != theGenericIon) {
        G4ProcessManager* pm =  theGenericIon->GetProcessManager();
        G4ProcessVector* v = pm->GetAlongStepProcessVector();
        size_t n = v->size();
        for(size_t j=0; j<n; ++j) {
          if((*v)[j] == this) {
            particle = theGenericIon;
            break;
          } 
        }
      }
    }
  }

  if( particle != &part ) {
    if(!isIon) {
      lManager->RegisterExtraParticle(&part, this);
    }
    if(1 < verboseLevel) {
      G4cout << "### G4VEnergyLossProcess::PreparePhysicsTable()"
             << " interrupted for "
             << part.GetParticleName() << "  isIon= " << isIon 
             << "  particle " << particle << "  GenericIon " << theGenericIon 
             << G4endl;
    }
    return;
  }

  Clean();
  lManager->PreparePhysicsTable(&part, this, isMaster);
  G4LossTableBuilder* bld = lManager->GetTableBuilder();

  // Base particle and set of models can be defined here
  InitialiseEnergyLossProcess(particle, baseParticle);

  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  size_t n = theCoupleTable->GetTableSize();

  theDEDXAtMaxEnergy.resize(n, 0.0);
  theRangeAtMaxEnergy.resize(n, 0.0);
  theEnergyOfCrossSectionMax.resize(n, 0.0);
  theCrossSectionMax.resize(n, DBL_MAX);

  // parameters of the process
  if(!actIntegral) { integral = theParameters->Integral(); }
  if(!actLossFluc) { lossFluctuationFlag = theParameters->LossFluctuation(); }
  rndmStepFlag = theParameters->UseCutAsFinalRange();
  if(!actMinKinEnergy) { minKinEnergy = theParameters->MinKinEnergy(); }
  if(!actMaxKinEnergy) { maxKinEnergy = theParameters->MaxKinEnergy(); }
  if(!actBinning) { 
    nBins = theParameters->NumberOfBinsPerDecade()
      *G4lrint(std::log10(maxKinEnergy/minKinEnergy));
  }
  maxKinEnergyCSDA = theParameters->MaxEnergyForCSDARange();
  nBinsCSDA = theParameters->NumberOfBinsPerDecade()
    *G4lrint(std::log10(maxKinEnergyCSDA/minKinEnergy));
  if(!actLinLossLimit) { linLossLimit = theParameters->LinearLossLimit(); }
  lambdaFactor = theParameters->LambdaFactor();
  if(isMaster) { SetVerboseLevel(theParameters->Verbose()); }
  else {  SetVerboseLevel(theParameters->WorkerVerbose()); }

  G4bool isElec = true;
  if(particle->GetPDGMass() > CLHEP::MeV) { isElec = false; }
  theParameters->DefineRegParamForLoss(this, isElec);

  G4double initialCharge = particle->GetPDGCharge();
  G4double initialMass   = particle->GetPDGMass();

  if (baseParticle) {
    massRatio = (baseParticle->GetPDGMass())/initialMass;
    G4double q = initialCharge/baseParticle->GetPDGCharge();
    chargeSqRatio = q*q;
    if(chargeSqRatio > 0.0) { reduceFactor = 1.0/(chargeSqRatio*massRatio); }
  }
  if(initialMass < MeV) {
    lowestKinEnergy = theParameters->LowestElectronEnergy();
  } else {
    lowestKinEnergy = theParameters->LowestMuHadEnergy();
  }

  // Tables preparation
  if (isMaster && !baseParticle) {

    if(theDEDXTable && isIonisation) {
      if(theIonisationTable && theDEDXTable != theIonisationTable) {
        theDEDXTable->clearAndDestroy();
        delete theDEDXTable;
        theDEDXTable = theIonisationTable;
      }   
      if(theDEDXSubTable && theIonisationSubTable && 
         theDEDXSubTable != theIonisationSubTable) {
        theDEDXSubTable->clearAndDestroy();
        delete theDEDXSubTable;
        theDEDXSubTable = theIonisationSubTable;
      }   
    }
    
    theDEDXTable = G4PhysicsTableHelper::PreparePhysicsTable(theDEDXTable);
    bld->InitialiseBaseMaterials(theDEDXTable);

    if(theDEDXSubTable) {
      theDEDXSubTable = 
        G4PhysicsTableHelper::PreparePhysicsTable(theDEDXSubTable);
    }

    if (theParameters->BuildCSDARange()) {
      theDEDXunRestrictedTable = 
        G4PhysicsTableHelper::PreparePhysicsTable(theDEDXunRestrictedTable);
      theCSDARangeTable = 
        G4PhysicsTableHelper::PreparePhysicsTable(theCSDARangeTable);
    }

    theLambdaTable = G4PhysicsTableHelper::PreparePhysicsTable(theLambdaTable);

    if(isIonisation) {
      theRangeTableForLoss = 
        G4PhysicsTableHelper::PreparePhysicsTable(theRangeTableForLoss);
      theInverseRangeTable = 
        G4PhysicsTableHelper::PreparePhysicsTable(theInverseRangeTable);  
    }

    if (nSCoffRegions && !lManager->SubCutProducer()) {
      theDEDXSubTable = 
        G4PhysicsTableHelper::PreparePhysicsTable(theDEDXSubTable);
      theSubLambdaTable = 
        G4PhysicsTableHelper::PreparePhysicsTable(theSubLambdaTable);
    }
  }
  /*
  G4cout << "** G4VEnergyLossProcess::PreparePhysicsTable() for " 
         << GetProcessName() << " and " << particle->GetParticleName()
         << " isMaster: " << isMaster << " isIonisation: " 
         << isIonisation << G4endl;
  G4cout << " theDEDX: " << theDEDXTable 
         << " theRange: " << theRangeTableForLoss
         << " theInverse: " << theInverseRangeTable
         << " theLambda: " << theLambdaTable << G4endl;
  */
  // forced biasing
  if(biasManager) { 
    biasManager->Initialise(part,GetProcessName(),verboseLevel); 
    biasFlag = false; 
  }

  // defined ID of secondary particles
  if(isMaster) {
    G4String nam1 = GetProcessName();
    G4String nam4 = nam1 + "_split";
    G4String nam5 = nam1 + "_subcut";
    secID   = G4PhysicsModelCatalog::Register(nam1); 
    biasID  = G4PhysicsModelCatalog::Register(nam4); 
    subsecID= G4PhysicsModelCatalog::Register(nam5);
  } 

  // initialisation of models
  G4int nmod = modelManager->NumberOfModels();
  for(G4int i=0; i<nmod; ++i) {
    G4VEmModel* mod = modelManager->GetModel(i);
    mod->SetMasterThread(isMaster);
    mod->SetAngularGeneratorFlag(
      theParameters->UseAngularGeneratorForIonisation());
    if(mod->HighEnergyLimit() > maxKinEnergy) {
      mod->SetHighEnergyLimit(maxKinEnergy);
    }
  }
  theCuts = modelManager->Initialise(particle, secondaryParticle, 
                                     theParameters->MinSubRange(), 
                                     verboseLevel);

  // Sub Cutoff 
  if(nSCoffRegions > 0) {
    if(theParameters->MinSubRange() < 1.0) { useSubCutoff = true; }

    theSubCuts = modelManager->SubCutoff();

    idxSCoffRegions = new G4bool[n]; 
    for (size_t j=0; j<n; ++j) {

      const G4MaterialCutsCouple* couple = 
        theCoupleTable->GetMaterialCutsCouple(j);
      const G4ProductionCuts* pcuts = couple->GetProductionCuts();
      
      G4bool reg = false;
      for(G4int i=0; i<nSCoffRegions; ++i) {
        if( pcuts == scoffRegions[i]->GetProductionCuts()) { 
          reg = true;
          break; 
        }
      }
      idxSCoffRegions[j] = reg;
    }
  }

  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::PrepearPhysicsTable() is done "
           << " for local " << particle->GetParticleName()
           << " isIon= " << isIon;
    if(baseParticle) { 
      G4cout << "; base: " << baseParticle->GetParticleName(); 
    }
    G4cout << " chargeSqRatio= " << chargeSqRatio
           << " massRatio= " << massRatio
           << " reduceFactor= " << reduceFactor << G4endl;
    if (nSCoffRegions) {
      G4cout << " SubCutoff Regime is ON for regions: " << G4endl;
      for (G4int i=0; i<nSCoffRegions; ++i) {
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
    if(baseParticle) { 
      G4cout << "; base: " << baseParticle->GetParticleName(); 
    }
    G4cout << " TablesAreBuilt= " << tablesAreBuilt
           << " isIon= " << isIon << "  " << this << G4endl;
  }

  if(&part == particle) {

    G4LossTableBuilder* bld = lManager->GetTableBuilder();
    if(isMaster) {
      theDensityFactor = bld->GetDensityFactors();
      theDensityIdx = bld->GetCoupleIndexes();
      lManager->BuildPhysicsTable(particle, this);

    } else {

      const G4VEnergyLossProcess* masterProcess = 
        static_cast<const G4VEnergyLossProcess*>(GetMasterProcess());

      // define density factors for worker thread
      bld->InitialiseBaseMaterials(masterProcess->DEDXTable()); 
      theDensityFactor = bld->GetDensityFactors();
      theDensityIdx = bld->GetCoupleIndexes();

      // copy table pointers from master thread
      SetDEDXTable(masterProcess->DEDXTable(),fRestricted);
      SetDEDXTable(masterProcess->DEDXTableForSubsec(),fSubRestricted);
      SetDEDXTable(masterProcess->DEDXunRestrictedTable(),fTotal);
      SetDEDXTable(masterProcess->IonisationTable(),fIsIonisation);
      SetDEDXTable(masterProcess->IonisationTableForSubsec(),fIsSubIonisation);
      SetRangeTableForLoss(masterProcess->RangeTableForLoss());
      SetCSDARangeTable(masterProcess->CSDARangeTable());
      SetSecondaryRangeTable(masterProcess->SecondaryRangeTable());
      SetInverseRangeTable(masterProcess->InverseRangeTable());
      SetLambdaTable(masterProcess->LambdaTable());
      SetSubLambdaTable(masterProcess->SubLambdaTable());
      isIonisation = masterProcess->IsIonisationProcess();

      tablesAreBuilt = true;  
      // local initialisation of models
      G4bool printing = true;
      G4int numberOfModels = modelManager->NumberOfModels();
      for(G4int i=0; i<numberOfModels; ++i) {
        G4VEmModel* mod = GetModelByIndex(i, printing);
        G4VEmModel* mod0= masterProcess->GetModelByIndex(i,printing);
        mod->InitialiseLocal(particle, mod0);
      }

      lManager->LocalPhysicsTables(particle, this);
    }
   
    // needs to be done only once
    safetyHelper->InitialiseHelper();
  }
  // explicitly defined printout by particle name
  G4String num = part.GetParticleName();
  if(1 < verboseLevel || 
     (0 < verboseLevel && (num == "e-" || 
                           num == "e+"    || num == "mu+" || 
                           num == "mu-"   || num == "proton"|| 
                           num == "pi+"   || num == "pi-" || 
                           num == "kaon+" || num == "kaon-" || 
                           num == "alpha" || num == "anti_proton" || 
                           num == "GenericIon")))
    { 
      PrintInfoDefinition(part); 
    }

  // Added tracking cut to avoid tracking artifacts
  // identify deexcitation flag
  if(isIonisation) { 
    //VI: seems not needed fParticleChange.SetLowEnergyLimit(lowestKinEnergy); 
    atomDeexcitation = lManager->AtomDeexcitation();
    if(nSCoffRegions > 0) { subcutProducer = lManager->SubCutProducer(); }
    if(atomDeexcitation) { 
      if(atomDeexcitation->IsPIXEActive()) { useDeexcitation = true; } 
    }
  }
  /*  
  G4cout << "** G4VEnergyLossProcess::BuildPhysicsTable() for " 
         << GetProcessName() << " and " << particle->GetParticleName()
         << " isMaster: " << isMaster << " isIonisation: " 
         << isIonisation << G4endl;
  G4cout << " theDEDX: " << theDEDXTable 
         << " theRange: " << theRangeTableForLoss
         << " theInverse: " << theInverseRangeTable
         << " theLambda: " << theLambdaTable << G4endl;
  */
  //if(1 < verboseLevel || verb) {
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
  if(1 < verboseLevel ) {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTable() of type " << tType
           << " for " << GetProcessName()
           << " and particle " << particle->GetParticleName()
           << G4endl;
  }
  G4PhysicsTable* table = nullptr;
  G4double emax = maxKinEnergy;
  G4int bin = nBins;

  if(fTotal == tType) {
    emax  = maxKinEnergyCSDA;
    bin   = nBinsCSDA;
    table = theDEDXunRestrictedTable;
  } else if(fRestricted == tType) {
    table = theDEDXTable;
  } else if(fSubRestricted == tType) {    
    table = theDEDXSubTable;
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
           << " table= " << table << "  " << this 
           << G4endl;
  }
  if(!table) { return table; }

  G4LossTableBuilder* bld = lManager->GetTableBuilder();
  G4bool splineFlag = theParameters->Spline();
  G4PhysicsLogVector* aVector = nullptr;
  G4PhysicsLogVector* bVector = nullptr;

  for(size_t i=0; i<numOfCouples; ++i) {

    if(1 < verboseLevel) {
      G4cout << "G4VEnergyLossProcess::BuildDEDXVector Idx= " << i 
             << "  flagTable=  " << table->GetFlag(i) 
             << " Flag= " << bld->GetFlag(i) << G4endl;
    }
    if(bld->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = 
        theCoupleTable->GetMaterialCutsCouple(i);
      if((*table)[i]) { delete (*table)[i]; }
      if(bVector) {
        aVector = new G4PhysicsLogVector(*bVector);
      } else {
        bVector = new G4PhysicsLogVector(minKinEnergy, emax, bin);
        aVector = bVector;
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
    if(2 < verboseLevel) G4cout << (*table) << G4endl;
  }

  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossProcess::BuildLambdaTable(G4EmTableType tType)
{
  G4PhysicsTable* table = nullptr;

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

  G4LossTableBuilder* bld = lManager->GetTableBuilder();
  theDensityFactor = bld->GetDensityFactors();
  theDensityIdx = bld->GetCoupleIndexes();

  G4bool splineFlag = theParameters->Spline();
  G4PhysicsLogVector* aVector = nullptr;
  G4double scale = G4Log(maxKinEnergy/minKinEnergy);

  for(size_t i=0; i<numOfCouples; ++i) {

    if (bld->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = 
        theCoupleTable->GetMaterialCutsCouple(i);
      delete (*table)[i];

      G4bool startNull = true;
      G4double emin = 
        MinPrimaryEnergy(particle,couple->GetMaterial(),(*theCuts)[i]);
      if(minKinEnergy > emin) { 
        emin = minKinEnergy; 
        startNull = false;
      }

      G4double emax = maxKinEnergy;
      if(emax <= emin) { emax = 2*emin; }
      G4int bin = G4lrint(nBins*G4Log(emax/emin)/scale);
      bin = std::max(bin, 3);
      aVector = new G4PhysicsLogVector(emin, emax, bin);
      aVector->SetSpline(splineFlag);

      modelManager->FillLambdaVector(aVector, couple, startNull, tType);
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

void 
G4VEnergyLossProcess::PrintInfoDefinition(const G4ParticleDefinition& part)
{
  if(0 < verboseLevel) {
    G4cout << std::setprecision(6);
    G4cout << G4endl << GetProcessName() << ":   for  "
           << part.GetParticleName()
           << "    SubType= " << GetProcessSubType() 
           << G4endl;
    G4cout << "      dE/dx and range tables from "
            << G4BestUnit(minKinEnergy,"Energy")
           << " to " << G4BestUnit(maxKinEnergy,"Energy")
           << " in " << nBins << " bins" << G4endl
           << "      Lambda tables from threshold to "
           << G4BestUnit(maxKinEnergy,"Energy")
           << ", " << theParameters->NumberOfBinsPerDecade() 
           << " bins per decade, spline: " 
           << theParameters->Spline()
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
      G4cout << "      SubLambdaTable address= " << theSubLambdaTable 
             << G4endl;
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
  if (!reg) {
    reg = regionStore->GetRegion("DefaultRegionForTheWorld", false);
  }

  // the region is in the list
  if (nSCoffRegions > 0) {
    for (G4int i=0; i<nSCoffRegions; ++i) {
      if (reg == scoffRegions[i]) {
        return;
      }
    }
  }
  // new region 
  if(val) {
    scoffRegions.push_back(reg);
    ++nSCoffRegions;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::StartTracking(G4Track* track)
{
  /*      
    G4cout << track->GetDefinition()->GetParticleName() 
           << " e(MeV)= " << track->GetKineticEnergy()
           << "  baseParticle " << baseParticle << " proc " << this;
    if(particle) G4cout << "  " << particle->GetParticleName();
    G4cout << " isIon= " << isIon << " dedx " << theDEDXTable <<G4endl;
  */
  // reset parameters for the new track
  theNumberOfInteractionLengthLeft = -1.0;
  currentInteractionLength = mfpKinEnergy = DBL_MAX; 
  preStepRangeEnergy = 0.0;

  // reset ion
  if(isIon) {
    chargeSqRatio = 0.5;

    G4double newmass = track->GetDefinition()->GetPDGMass();
    if(baseParticle) {
      massRatio = baseParticle->GetPDGMass()/newmass;
    } else if(theGenericIon) {
      massRatio = proton_mass_c2/newmass;
    } else {
      massRatio = 1.0;
    }
  }  
  // forced biasing only for primary particles
  if(biasManager) {
    if(0 == track->GetParentID()) {
      // primary particle
      biasFlag = true; 
      biasManager->ResetForcedInteraction(); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::AlongStepGetPhysicalInteractionLength(
                             const G4Track&,G4double,G4double,G4double&,
                             G4GPILSelection* selection)
{
  G4double x = DBL_MAX;
  *selection = aGPILSelection;
  if(isIonisation && currentModel->IsActive(preStepScaledEnergy)) {
    fRange = GetScaledRangeForScaledEnergy(preStepScaledEnergy)*reduceFactor;
    x = fRange;
    G4double finR = finalRange;
    if(rndmStepFlag) { 
      finR = std::min(finR,
                      currentCouple->GetProductionCuts()->GetProductionCut(1));
    }
    if(fRange > finR) { 
      x = fRange*dRoverRange + finR*(1.0 - dRoverRange)*(2.0 - finR/fRange); 
    }
    
   // if(particle->GetPDGMass() > 0.9*GeV)
    /*
    G4cout<<GetProcessName()<<": e= "<<preStepKinEnergy
          <<" range= "<<fRange << " idx= " << basedCoupleIndex
              << " finR= " << finR
          << " limit= " << x <<G4endl;
    G4cout << "massRatio= " << massRatio << " Q^2= " << chargeSqRatio 
           << " finR= " << finR << " dRoverRange= " << dRoverRange 
           << " finalRange= " << finalRange << G4endl;
    */
  }
  //G4cout<<GetProcessName()<<": e= "<<preStepKinEnergy 
  //<<" stepLimit= "<<x<<G4endl;
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

  // initialisation of material, mass, charge, model 
  // at the beginning of the step
  DefineMaterial(track.GetMaterialCutsCouple());
  preStepKinEnergy    = track.GetKineticEnergy();
  preStepScaledEnergy = preStepKinEnergy*massRatio;
  SelectModel(preStepScaledEnergy);

  if(!currentModel->IsActive(preStepScaledEnergy)) { 
    theNumberOfInteractionLengthLeft = -1.0;
    currentInteractionLength = DBL_MAX;
    return x; 
  }

  // change effective charge of an ion on fly
  if(isIon) {
    G4double q2 = currentModel->ChargeSquareRatio(track);
    if(q2 != chargeSqRatio && q2 > 0.0) {
      chargeSqRatio = q2;
      fFactor = q2*biasFactor*(*theDensityFactor)[currentCoupleIndex];
      reduceFactor = 1.0/(fFactor*massRatio);
    }
  }
  //  if(particle->GetPDGMass() > 0.9*GeV)
  //G4cout << "q2= "<<chargeSqRatio << " massRatio= " << massRatio << G4endl; 

  // forced biasing only for primary particles
  if(biasManager) {
    if(0 == track.GetParentID()) {
      if(biasFlag && 
         biasManager->ForcedInteractionRegion(currentCoupleIndex)) {
        return biasManager->GetStepLimit(currentCoupleIndex, previousStepSize);
      }
    }
  }

  // compute mean free path
  if(preStepScaledEnergy < mfpKinEnergy) {
    if (integral) { ComputeLambdaForScaledEnergy(preStepScaledEnergy); }
    else  { preStepLambda = GetLambdaForScaledEnergy(preStepScaledEnergy); }

    // zero cross section
    if(preStepLambda <= 0.0) { 
      theNumberOfInteractionLengthLeft = -1.0;
      currentInteractionLength = DBL_MAX;
    }
  }

  // non-zero cross section
  if(preStepLambda > 0.0) { 
    if (theNumberOfInteractionLengthLeft < 0.0) {

      // beggining of tracking (or just after DoIt of this process)
      theNumberOfInteractionLengthLeft =  -G4Log( G4UniformRand() );
      theInitialNumberOfInteractionLength = theNumberOfInteractionLengthLeft; 

    } else if(currentInteractionLength < DBL_MAX) {

      // subtract NumberOfInteractionLengthLeft using previous step
      theNumberOfInteractionLengthLeft -= 
        previousStepSize/currentInteractionLength;

      theNumberOfInteractionLengthLeft = 
        std::max(theNumberOfInteractionLengthLeft, 0.0);
    }

    // new mean free path and step limit
    currentInteractionLength = 1.0/preStepLambda;
    x = theNumberOfInteractionLengthLeft * currentInteractionLength;
  }
#ifdef G4VERBOSE
  if (verboseLevel>2){
    //  if(particle->GetPDGMass() > 0.9*GeV){
    G4cout << "G4VEnergyLossProcess::PostStepGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" << G4endl; 
    G4cout << " for " << track.GetDefinition()->GetParticleName() 
           << " in Material  " <<  currentMaterial->GetName()
           << " Ekin(MeV)= " << preStepKinEnergy/MeV 
           << "  " << track.GetMaterial()->GetName()
           <<G4endl;
    G4cout << "MeanFreePath = " << currentInteractionLength/cm << "[cm]" 
           << "InteractionLength= " << x/cm <<"[cm] " <<G4endl;
  }
#endif
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
  if(length <= 0.0) { return &fParticleChange; }
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
           << "  rf= " << reduceFactor
           << "  q^2= " << chargeSqRatio
           << " md= " << d->GetPDGMass()
           << "  status= " << track.GetTrackStatus()
           << "  " << track.GetMaterial()->GetName()
           << G4endl;
  }
  */

  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();

  // define new weight for primary and secondaries
  G4double weight = fParticleChange.GetParentWeight();
  if(weightFlag) {
    weight /= biasFactor;
    fParticleChange.ProposeWeight(weight);
  }

  // stopping
  if (length >= fRange || preStepKinEnergy <= lowestKinEnergy) {
    eloss = preStepKinEnergy;
    if (useDeexcitation) {
      atomDeexcitation->AlongStepDeexcitation(scTracks, step, 
                                              eloss, currentCoupleIndex);
      if(scTracks.size() > 0) { FillSecondariesAlongStep(eloss, weight); }
      eloss = std::max(eloss, 0.0);
    }
    fParticleChange.SetProposedKineticEnergy(0.0);
    fParticleChange.ProposeLocalEnergyDeposit(eloss);
    return &fParticleChange;
  }
  //G4cout << theDEDXTable << "  idx= " << basedCoupleIndex 
  // << "  " << GetProcessName() << "  "<< currentMaterial->GetName()<<G4endl;
  //if(particle->GetParticleName() == "e-")G4cout << (*theDEDXTable) <<G4endl;
  // Short step
  eloss = GetDEDXForScaledEnergy(preStepScaledEnergy)*length;

  //G4cout << "eloss= " << eloss << G4endl;

  // Long step
  if(eloss > preStepKinEnergy*linLossLimit) {

    G4double x = (fRange - length)/reduceFactor;
    //G4cout << "x= " << x << "  " << theInverseRangeTable << G4endl;
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

  G4double cut  = (*theCuts)[currentCoupleIndex];
  G4double esec = 0.0;

  //G4cout << "cut= " << cut << " useSubCut= " << useSubCutoff << G4endl;

  // SubCutOff 
  if(useSubCutoff && !subcutProducer) {
    if(idxSCoffRegions[currentCoupleIndex]) {

      G4bool yes = false;
      G4StepPoint* prePoint = step.GetPreStepPoint();

      // Check boundary
      if(prePoint->GetStepStatus() == fGeomBoundary) { yes = true; }

      // Check PrePoint
      else {
        G4double preSafety  = prePoint->GetSafety();
        G4double rcut = 
          currentCouple->GetProductionCuts()->GetProductionCut(1);

        // recompute presafety
        if(preSafety < rcut) {
          preSafety = safetyHelper->ComputeSafety(prePoint->GetPosition(),
                                                  rcut);
        }

        if(preSafety < rcut) { yes = true; }

        // Check PostPoint
        else {
          G4double postSafety = preSafety - length; 
          if(postSafety < rcut) {
            postSafety = safetyHelper->ComputeSafety(
              step.GetPostStepPoint()->GetPosition(), rcut);
            if(postSafety < rcut) { yes = true; }
          }
        }
      }
  
      // Decided to start subcut sampling
      if(yes) {

        cut = (*theSubCuts)[currentCoupleIndex];
         eloss -= GetSubDEDXForScaledEnergy(preStepScaledEnergy)*length;
        esec = SampleSubCutSecondaries(scTracks, step, 
                                       currentModel,currentCoupleIndex);
        // add bremsstrahlung sampling
        /*
        if(nProcesses > 0) {
          for(G4int i=0; i<nProcesses; ++i) {
            (scProcesses[i])->SampleSubCutSecondaries(
                scTracks, step, (scProcesses[i])->
                SelectModelForMaterial(preStepKinEnergy, currentCoupleIndex),
                currentCoupleIndex);
          }
        } 
        */
      }   
    }
  }

  // Corrections, which cannot be tabulated
  if(isIon) {
    G4double eadd = 0.0;
    G4double eloss_before = eloss;
    currentModel->CorrectionsAlongStep(currentCouple, dynParticle, 
                                       eloss, eadd, length);
    if(eloss < 0.0) { eloss = 0.5*eloss_before; }
  }

  // Sample fluctuations
  if (lossFluctuationFlag) {
    G4VEmFluctuationModel* fluc = currentModel->GetModelOfFluctuations();
    if(eloss + esec < preStepKinEnergy) {

      G4double tmax = 
        std::min(currentModel->MaxSecondaryKinEnergy(dynParticle),cut);
      eloss = fluc->SampleFluctuations(currentCouple,dynParticle,
                                       tmax,length,eloss);
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
    G4double esecfluo = preStepKinEnergy - esec;
    G4double de = esecfluo;
    //G4double eloss0 = eloss;
    /*
    G4cout << "### 1: E(keV)= " << preStepKinEnergy/keV
           << " Efluomax(keV)= " << de/keV
           << " Eloss(keV)= " << eloss/keV << G4endl; 
    */
    atomDeexcitation->AlongStepDeexcitation(scTracks, step, 
                                            de, currentCoupleIndex);

    // sum of de-excitation energies
    esecfluo -= de;

    // subtracted from energy loss
    if(eloss >= esecfluo) {
      esec  += esecfluo;
      eloss -= esecfluo;
    } else {
      esec += esecfluo;
      eloss = 0.0; 
    } 
    /*    
    if(esecfluo > 0.0) {
      G4cout << "### 2: E(keV)= " << preStepKinEnergy/keV
             << " Esec(keV)= " << esec/keV
             << " Esecf(kV)= " << esecfluo/keV
             << " Eloss0(kV)= " << eloss0/keV
             << " Eloss(keV)= " << eloss/keV 
             << G4endl; 
    } 
    */   
  }
  if(subcutProducer && idxSCoffRegions[currentCoupleIndex]) {
    subcutProducer->SampleSecondaries(step, scTracks, eloss, cut);
  }
  if(scTracks.size() > 0) { FillSecondariesAlongStep(eloss, weight); }

  // Energy balance
  G4double finalT = preStepKinEnergy - eloss - esec;
  if (finalT <= lowestKinEnergy) {
    eloss += finalT;
    finalT = 0.0;
  } else if(isIon) {
    fParticleChange.SetProposedCharge(
      currentModel->GetParticleCharge(track.GetParticleDefinition(),
                                      currentMaterial,finalT));
  }

  eloss = std::max(eloss, 0.0);

  fParticleChange.SetProposedKineticEnergy(finalT);
  fParticleChange.ProposeLocalEnergyDeposit(eloss);
  /*
  if(-1 < verboseLevel) {
    G4double del = finalT + eloss + esec - preStepKinEnergy;
    G4cout << "Final value eloss(MeV)= " << eloss/MeV
           << " preStepKinEnergy= " << preStepKinEnergy
           << " postStepKinEnergy= " << finalT
           << " de(keV)= " << del/keV
           << " lossFlag= " << lossFluctuationFlag
           << "  status= " << track.GetTrackStatus()
           << G4endl;
  }
  */
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4VEnergyLossProcess::FillSecondariesAlongStep(G4double&, G4double& weight)
{
  G4int n0 = scTracks.size();

  // weight may be changed by biasing manager
  if(biasManager) {
    if(biasManager->SecondaryBiasingRegion(currentCoupleIndex)) {
      weight *=
        biasManager->ApplySecondaryBiasing(scTracks, currentCoupleIndex);
    }
  } 

  // fill secondaries
  G4int n = scTracks.size();
  fParticleChange.SetNumberOfSecondaries(n);

  for(G4int i=0; i<n; ++i) {
    G4Track* t = scTracks[i];
    if(t) {
      t->SetWeight(weight); 
      pParticleChange->AddSecondary(t);
      if(i >= n0) { t->SetCreatorModelIndex(biasID); }
      //G4cout << "Secondary(along step) has weight " << t->GetWeight() 
      //<< ", kenergy " << t->GetKineticEnergy()/MeV << " MeV" <<G4endl;
    }
  }
  scTracks.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4VEnergyLossProcess::SampleSubCutSecondaries(std::vector<G4Track*>& tracks, 
                                              const G4Step& step, 
                                              G4VEmModel* model,
                                              G4int idx) 
{
  // Fast check weather subcutoff can work
  G4double esec = 0.0;
  G4double subcut = (*theSubCuts)[idx];
  G4double cut = (*theCuts)[idx];
  if(cut <= subcut) { return esec; }

  const G4Track* track = step.GetTrack();
  const G4DynamicParticle* dp = track->GetDynamicParticle();
  G4double e = dp->GetKineticEnergy()*massRatio;
  G4double cross = (*theDensityFactor)[idx]*chargeSqRatio
    *(((*theSubLambdaTable)[(*theDensityIdx)[idx]])->Value(e, idxSubLambda));
  G4double length = step.GetStepLength();

  // negligible probability to get any interaction
  if(length*cross < perMillion) { return esec; }
  /*      
  if(-1 < verboseLevel) 
    G4cout << "<<< Subcutoff for " << GetProcessName()
           << " cross(1/mm)= " << cross*mm << ">>>"
           << " e(MeV)= " << preStepScaledEnergy
           << " matIdx= " << currentCoupleIndex
           << G4endl;
  */

  // Sample subcutoff secondaries
  G4StepPoint* preStepPoint = step.GetPreStepPoint();
  G4StepPoint* postStepPoint = step.GetPostStepPoint();
  G4ThreeVector prepoint = preStepPoint->GetPosition();
  G4ThreeVector dr = postStepPoint->GetPosition() - prepoint;
  G4double pretime = preStepPoint->GetGlobalTime();
  G4double dt = postStepPoint->GetGlobalTime() - pretime;
  G4double fragment = 0.0;

  do {
    G4double del = -G4Log(G4UniformRand())/cross;
    fragment += del/length;
    if (fragment > 1.0) { break; }

    // sample secondaries
    secParticles.clear();
    model->SampleSecondaries(&secParticles,track->GetMaterialCutsCouple(),
                             dp,subcut,cut);

    // position of subcutoff particles
    G4ThreeVector r = prepoint + fragment*dr;
    std::vector<G4DynamicParticle*>::iterator it;
    for(it=secParticles.begin(); it!=secParticles.end(); ++it) {

      G4Track* t = new G4Track((*it), pretime + fragment*dt, r);
      t->SetTouchableHandle(track->GetTouchableHandle());
      t->SetCreatorModelIndex(subsecID);
      tracks.push_back(t);
      esec += t->GetKineticEnergy();
      if (t->GetParticleDefinition() == thePositron) { 
        esec += 2.0*electron_mass_c2; 
      }

        /*        
        if(-1 < verboseLevel) 
          G4cout << "New track " 
                 << t->GetParticleDefinition()->GetParticleName()
                 << " e(keV)= " << t->GetKineticEnergy()/keV
                 << " fragment= " << fragment
                 << G4endl;
        */
    }
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while (fragment <= 1.0);
  return esec;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEnergyLossProcess::PostStepDoIt(const G4Track& track,
                                                      const G4Step& step)
{
  // In all cases clear number of interaction lengths
  theNumberOfInteractionLengthLeft = -1.0;
  mfpKinEnergy = currentInteractionLength = DBL_MAX; 

  fParticleChange.InitializeForPostStep(track);
  G4double finalT = track.GetKineticEnergy();
  if(finalT <= lowestKinEnergy) { return &fParticleChange; }

  G4double postStepScaledEnergy = finalT*massRatio;

  if(!currentModel->IsActive(postStepScaledEnergy)) { 
    return &fParticleChange; 
  }
  /*
  if(-1 < verboseLevel) {
    G4cout << GetProcessName()
           << "::PostStepDoIt: E(MeV)= " << finalT/MeV
           << G4endl;
  }
  */

  // forced process - should happen only once per track
  if(biasFlag) {
    if(biasManager->ForcedInteractionRegion(currentCoupleIndex)) {
      biasFlag = false;
    }
  }

  // Integral approach
  if (integral) {
    G4double lx = GetLambdaForScaledEnergy(postStepScaledEnergy);
    /*
    if(preStepLambda<lx && 1 < verboseLevel) {
      G4cout << "WARNING: for " << particle->GetParticleName()
             << " and " << GetProcessName()
             << " E(MeV)= " << finalT/MeV
             << " preLambda= " << preStepLambda 
             << " < " << lx << " (postLambda) "
             << G4endl;
    }
    */
    if(lx <= 0.0 || preStepLambda*G4UniformRand() > lx) {
      return &fParticleChange;
    }
  }

  SelectModel(postStepScaledEnergy);

  // define new weight for primary and secondaries
  G4double weight = fParticleChange.GetParentWeight();
  if(weightFlag) {
    weight /= biasFactor;
    fParticleChange.ProposeWeight(weight);
  }

  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  G4double tcut = (*theCuts)[currentCoupleIndex];

  // sample secondaries
  secParticles.clear();
  //G4cout<< "@@@ Eprimary= "<<dynParticle->GetKineticEnergy()/MeV
  //        << " cut= " << tcut/MeV << G4endl;
  currentModel->SampleSecondaries(&secParticles, currentCouple, 
                                  dynParticle, tcut);

  G4int num0 = secParticles.size();

  // bremsstrahlung splitting or Russian roulette  
  if(biasManager) {
    if(biasManager->SecondaryBiasingRegion(currentCoupleIndex)) {
      G4double eloss = 0.0;
      weight *= biasManager->ApplySecondaryBiasing(
                                      secParticles,
                                      track, currentModel, 
                                      &fParticleChange, eloss,
                                      currentCoupleIndex, tcut, 
                                      step.GetPostStepPoint()->GetSafety());
      if(eloss > 0.0) {
        eloss += fParticleChange.GetLocalEnergyDeposit();
        fParticleChange.ProposeLocalEnergyDeposit(eloss);
      }
    }
  }

  // save secondaries
  G4int num = secParticles.size();
  if(num > 0) {

    fParticleChange.SetNumberOfSecondaries(num);
    G4double time = track.GetGlobalTime();

    for (G4int i=0; i<num; ++i) {
      if(secParticles[i]) {
        G4Track* t = new G4Track(secParticles[i], time, track.GetPosition());
        t->SetTouchableHandle(track.GetTouchableHandle());
        t->SetWeight(weight); 
        if(i < num0) { t->SetCreatorModelIndex(secID); }
        else         { t->SetCreatorModelIndex(biasID); }

        //G4cout << "Secondary(post step) has weight " << t->GetWeight() 
        //       << ", kenergy " << t->GetKineticEnergy()/MeV << " MeV" 
        //       << " time= " << time/ns << " ns " << G4endl;
        pParticleChange->AddSecondary(t);
      }
    }
  }

  if(0.0 == fParticleChange.GetProposedKineticEnergy() &&
     fAlive == fParticleChange.GetTrackStatus()) {
    if(particle->GetProcessManager()->GetAtRestProcessVector()->size() > 0)
         { fParticleChange.ProposeTrackStatus(fStopButAlive); }
    else { fParticleChange.ProposeTrackStatus(fStopAndKill); }
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
  //G4cout << "G4VEnergyLossProcess::StorePhysicsTable: " << part->GetParticleName()
  //         << "  " << directory << "  " << ascii << G4endl;
  if (!isMaster || baseParticle || part != particle ) return res;

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

  if ( !res ) {
    if(1 < verboseLevel) {
      G4cout << "Physics tables are stored for " 
             << particle->GetParticleName()
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
  if (!isMaster) return res;
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

      if(!RetrieveTable(part,theDEDXunRestrictedTable,ascii,directory,
                        "DEDXnr",false)) 
        {res = false;}

      if(!RetrieveTable(part,theCSDARangeTable,ascii,directory,
                        "CSDARange",false)) 
        {res = false;}

      if(!RetrieveTable(part,theInverseRangeTable,ascii,directory,
                        "InverseRange",fpi)) 
        {res = false;}

      if(!RetrieveTable(part,theLambdaTable,ascii,directory,"Lambda",true)) 
        {res = false;}

      G4bool yes = false;
      if(nSCoffRegions > 0) {yes = true;}

      if(!RetrieveTable(part,theDEDXSubTable,ascii,directory,"SubDEDX",yes)) 
        {res = false;}

      if(!RetrieveTable(part,theSubLambdaTable,ascii,directory,
                        "SubLambda",yes)) 
        {res = false;}

      if(!fpi) yes = false;
      if(!RetrieveTable(part,theIonisationSubTable,ascii,directory,
                        "SubIonisation",yes))
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
  //G4cout << "G4VEnergyLossProcess::StoreTable: " << aTable
  //         << "  " << directory << "  " << tname << G4endl;
  G4bool res = true;
  if ( aTable ) {
    const G4String name = GetPhysicsTableFileName(part,directory,tname,ascii);
    G4cout << name << G4endl;
    //G4cout << *aTable << G4endl;
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
        if(theParameters->Spline()) {
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
  tmax = std::min(tmax,(*theCuts)[currentCoupleIndex]);
  G4double d = 0.0;
  G4VEmFluctuationModel* fm = currentModel->GetModelOfFluctuations();
  if(fm) { d = fm->Dispersion(currentMaterial,dp,tmax,length); }
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
    cross = GetLambdaForScaledEnergy(kineticEnergy*massRatio);
  } else {
    SelectModel(kineticEnergy*massRatio);
    cross = biasFactor*(*theDensityFactor)[currentCoupleIndex]
      *(currentModel->CrossSectionPerVolume(currentMaterial,
                                            particle, kineticEnergy,
                                            (*theCuts)[currentCoupleIndex]));
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
  if(0.0 < preStepLambda) { x = 1.0/preStepLambda; }
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
G4VEnergyLossProcess::LambdaPhysicsVector(const G4MaterialCutsCouple*, 
                                          G4double)
{
  G4PhysicsVector* v = 
    new G4PhysicsLogVector(minKinEnergy, maxKinEnergy, nBins);
  v->SetSpline(theParameters->Spline());
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

void 
G4VEnergyLossProcess::SetDEDXTable(G4PhysicsTable* p, G4EmTableType tType)
{
  if(fTotal == tType) {
    theDEDXunRestrictedTable = p;
    if(p) {
      size_t n = p->length();
      G4PhysicsVector* pv = (*p)[0];
      G4double emax = maxKinEnergyCSDA;

      G4LossTableBuilder* bld = lManager->GetTableBuilder();
      theDensityFactor = bld->GetDensityFactors();
      theDensityIdx = bld->GetCoupleIndexes();

      for (size_t i=0; i<n; ++i) {
        G4double dedx = 0.0; 
        pv = (*p)[i];
        if(pv) { 
          dedx = pv->Value(emax, idxDEDXunRestricted); 
        } else {
          pv = (*p)[(*theDensityIdx)[i]];
          if(pv) { 
            dedx = 
              pv->Value(emax, idxDEDXunRestricted)*(*theDensityFactor)[i]; 
          }
        }
        theDEDXAtMaxEnergy[i] = dedx;
        //G4cout << "i= " << i << " emax(MeV)= " << emax/MeV<< " dedx= " 
        //     << dedx << G4endl;
      }
    }

  } else if(fRestricted == tType) {
    /*
      G4cout<< "G4VEnergyLossProcess::SetDEDXTable "
            << particle->GetParticleName()
            << " oldTable " << theDEDXTable << " newTable " << p 
            << " ion " << theIonisationTable 
            << " IsMaster " << isMaster 
            << " " << GetProcessName() << G4endl;
      G4cout << (*p) << G4endl;
    */
    theDEDXTable = p;
  } else if(fSubRestricted == tType) {
      theDEDXSubTable = p;
  } else if(fIsIonisation == tType) {
    /*
      G4cout<< "G4VEnergyLossProcess::SetIonisationTable "
            << particle->GetParticleName()
            << " oldTable " << theDEDXTable << " newTable " << p 
            << " ion " << theIonisationTable 
            << " IsMaster " << isMaster 
            << " " << GetProcessName() << G4endl;
    */
    theIonisationTable = p;
  } else if(fIsSubIonisation == tType) {
    theIonisationSubTable = p;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetCSDARangeTable(G4PhysicsTable* p)
{
  theCSDARangeTable = p; 

  if(p) {
    size_t n = p->length();
    G4PhysicsVector* pv;
    G4double emax = maxKinEnergyCSDA;

    for (size_t i=0; i<n; ++i) {
      pv = (*p)[i];
      G4double rmax = 0.0;
      if(pv) { rmax = pv->Value(emax, idxCSDA); }
      else {
        pv = (*p)[(*theDensityIdx)[i]];
        if(pv) { rmax = pv->Value(emax, idxCSDA)/(*theDensityFactor)[i]; }
      }
      theRangeAtMaxEnergy[i] = rmax;
      //G4cout << "i= " << i << " Emax(MeV)= " << emax/MeV << " Rmax= " 
      //<< rmax<< G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetRangeTableForLoss(G4PhysicsTable* p)
{
  theRangeTableForLoss = p;
  if(1 < verboseLevel) {
    G4cout << "### Set Range table " << p 
           << " for " << particle->GetParticleName()
           << " and process " << GetProcessName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetSecondaryRangeTable(G4PhysicsTable* p)
{
  theSecondaryRangeTable = p;
  if(1 < verboseLevel) {
    G4cout << "### Set SecondaryRange table " << p 
           << " for " << particle->GetParticleName()
           << " and process " << GetProcessName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetInverseRangeTable(G4PhysicsTable* p)
{
  theInverseRangeTable = p;
  if(1 < verboseLevel) {
    G4cout << "### Set InverseRange table " << p 
           << " for " << particle->GetParticleName()
           << " and process " << GetProcessName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetLambdaTable(G4PhysicsTable* p)
{
  if(1 < verboseLevel) {
    G4cout << "### Set Lambda table " << p 
           << " for " << particle->GetParticleName()
           << " and process " << GetProcessName() << G4endl;
    //G4cout << *p << G4endl;
  }
  theLambdaTable = p; 
  tablesAreBuilt = true;

  G4LossTableBuilder* bld = lManager->GetTableBuilder();
  theDensityFactor = bld->GetDensityFactors();
  theDensityIdx = bld->GetCoupleIndexes();

  if(theLambdaTable) {
    size_t n = theLambdaTable->length();
    G4PhysicsVector* pv = (*theLambdaTable)[0];
    G4double e, ss, smax, emax;

    size_t i;

    // first loop on existing vectors
    for (i=0; i<n; ++i) {
      pv = (*theLambdaTable)[i];
      if(pv) {
        size_t nb = pv->GetVectorLength();
        emax = DBL_MAX;
        smax = 0.0;
        if(nb > 0) {
          for (size_t j=0; j<nb; ++j) {
            e = pv->Energy(j);
            ss = (*pv)(j);
            if(ss > smax) {
              smax = ss;
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
    // second loop using base materials
    for (i=0; i<n; ++i) {
      pv = (*theLambdaTable)[i];
      if(!pv){
        G4int j = (*theDensityIdx)[i];
        theEnergyOfCrossSectionMax[i] = theEnergyOfCrossSectionMax[j];
        theCrossSectionMax[i] = (*theDensityFactor)[i]*theCrossSectionMax[j];
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetSubLambdaTable(G4PhysicsTable* p)
{
  theSubLambdaTable = p;
  if(1 < verboseLevel) {
    G4cout << "### Set SebLambda table " << p 
           << " for " << particle->GetParticleName()
           << " and process " << GetProcessName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Element* G4VEnergyLossProcess::GetCurrentElement() const
{
  const G4Element* elm = nullptr;
  if(currentModel) { elm = currentModel->GetCurrentElement(); }
  return elm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetCrossSectionBiasingFactor(G4double f, 
                                                        G4bool flag)
{
  if(f > 0.0) { 
    biasFactor = f; 
    weightFlag = flag;
    if(1 < verboseLevel) {
      G4cout << "### SetCrossSectionBiasingFactor: for " 
             << " process " << GetProcessName()
             << " biasFactor= " << f << " weightFlag= " << flag 
             << G4endl; 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4VEnergyLossProcess::ActivateForcedInteraction(G4double length, 
                                                const G4String& region,
                                                G4bool flag)
{
  if(!biasManager) { biasManager = new G4EmBiasingManager(); }
  if(1 < verboseLevel) {
    G4cout << "### ActivateForcedInteraction: for " 
           << " process " << GetProcessName()
           << " length(mm)= " << length/mm
           << " in G4Region <" << region
           << "> weightFlag= " << flag 
           << G4endl; 
  }
  weightFlag = flag;
  biasManager->ActivateForcedInteraction(length, region);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4VEnergyLossProcess::ActivateSecondaryBiasing(const G4String& region, 
                                               G4double factor, 
                                               G4double energyLimit)
{
  if (0.0 <= factor) {

    // Range cut can be applied only for e-
    if(0.0 == factor && secondaryParticle != G4Electron::Electron())
      { return; }

    if(!biasManager) { biasManager = new G4EmBiasingManager(); }
    biasManager->ActivateSecondaryBiasing(region, factor, energyLimit);
    if(1 < verboseLevel) {
      G4cout << "### ActivateSecondaryBiasing: for " 
             << " process " << GetProcessName()
             << " factor= " << factor
             << " in G4Region <" << region 
             << "> energyLimit(MeV)= " << energyLimit/MeV
             << G4endl; 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetIonisation(G4bool val)
{
  isIonisation = val;
  if(val) { aGPILSelection = CandidateForSelection; }
  else    { aGPILSelection = NotCandidateForSelection; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 void G4VEnergyLossProcess::SetLinearLossLimit(G4double val)
{
  if(0.0 < val && val < 1.0) { 
    linLossLimit = val;
    actLinLossLimit = true; 
  } else { PrintWarning("SetLinearLossLimit", val); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4VEnergyLossProcess::SetStepFunction(G4double v1, G4double v2, G4bool lock)
{
  if(actStepFunc) { return; }
  actStepFunc = lock;
  if(0.0 < v1 && 0.0 < v2 && v2 < 1.e+50) { 
    dRoverRange = std::min(1.0, v1);
    finalRange = v2;
  } else if(v1 <= 0.0) {
    PrintWarning("SetStepFunction", v1); 
  } else {
    PrintWarning("SetStepFunction", v2); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetLowestEnergyLimit(G4double val)
{
  if(1.e-18 < val && val < 1.e+50) { lowestKinEnergy = val; }
  else { PrintWarning("SetLowestEnergyLimit", val); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetDEDXBinning(G4int n)
{
  if(2 < n && n < 1000000000) { 
    nBins = n; 
    actBinning = true;
  } else {
    G4double e = (G4double)n;
    PrintWarning("SetDEDXBinning", e); 
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetMinKinEnergy(G4double e)
{
  if(1.e-18 < e && e < maxKinEnergy) { 
    minKinEnergy = e; 
    actMinKinEnergy = true;
  } else { PrintWarning("SetMinKinEnergy", e); } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetMaxKinEnergy(G4double e)
{
  if(minKinEnergy < e && e < 1.e+50) { 
    maxKinEnergy = e;
    actMaxKinEnergy = true;
    if(e < maxKinEnergyCSDA) { maxKinEnergyCSDA = e; }
  } else { PrintWarning("SetMaxKinEnergy", e); } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::PrintWarning(G4String tit, G4double val)
{
  G4String ss = "G4VEnergyLossProcess::" + tit; 
  G4ExceptionDescription ed;
  ed << "Parameter is out of range: " << val 
     << " it will have no effect!\n" << "  Process " 
     << GetProcessName() << "  nbins= " << nBins 
     << " Emin(keV)= " << minKinEnergy/keV 
     << " Emax(GeV)= " << maxKinEnergy/GeV;
  G4Exception(ss, "em0044", JustWarning, ed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::ProcessDescription(std::ostream& outFile) const
{
  outFile << "Energy loss process <" << GetProcessName() << ">" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
