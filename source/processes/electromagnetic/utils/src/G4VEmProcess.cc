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
// GEANT4 Class file
//
//
// File name:     G4VEmProcess
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 01.10.2003
//
// Modifications: by V.Ivanchenko
//
// Class Description: based class for discrete and rest/discrete EM processes
//

// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VEmProcess.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4LossTableManager.hh"
#include "G4LossTableBuilder.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEmModel.hh"
#include "G4DataVector.hh"
#include "G4PhysicsTable.hh"
#include "G4EmDataHandler.hh"
#include "G4PhysicsLogVector.hh"
#include "G4VParticleChange.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4EmBiasingManager.hh"
#include "G4GenericIon.hh"
#include "G4Log.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmProcess::G4VEmProcess(const G4String& name, G4ProcessType type):
  G4VDiscreteProcess(name, type),
  secondaryParticle(nullptr),
  buildLambdaTable(true),
  numberOfModels(0),
  theLambdaTable(nullptr),
  theLambdaTablePrim(nullptr),
  integral(false),
  applyCuts(false),
  startFromNull(false),
  splineFlag(true),
  isIon(false),
  currentCouple(nullptr),
  isTheMaster(true),
  masterProc(nullptr),
  theData(nullptr),
  currentModel(nullptr),
  particle(nullptr),
  currentParticle(nullptr)
{
  theParameters = G4EmParameters::Instance();
  SetVerboseLevel(1);

  // Size of tables assuming spline
  minKinEnergy = 0.1*keV;
  maxKinEnergy = 100.0*TeV;
  nLambdaBins  = 84;
  minKinEnergyPrim = DBL_MAX;
  actBinning = actSpline = actMinKinEnergy = actMaxKinEnergy = false;

  // default lambda factor
  lambdaFactor    = 0.8;
  logLambdaFactor = G4Log(lambdaFactor);

  // default limit on polar angle
  biasFactor = fFactor = 1.0;

  // particle types
  theGamma     = G4Gamma::Gamma();
  theElectron  = G4Electron::Electron();
  thePositron  = G4Positron::Positron();

  theCuts = theCutsGamma = theCutsElectron = theCutsPositron = nullptr;

  pParticleChange = &fParticleChange;
  fParticleChange.SetSecondaryWeightByProcess(true);
  secParticles.reserve(5);

  baseMaterial = currentMaterial = nullptr;

  preStepLambda = preStepKinEnergy = 0.0;
  preStepLogKinEnergy = LOG_EKIN_MIN;
  mfpKinEnergy  = DBL_MAX;
  massRatio     = 1.0;

  currentCoupleIndex = basedCoupleIndex = 0;

  modelManager = new G4EmModelManager();
  biasManager  = nullptr;
  biasFlag     = false; 
  weightFlag   = false;
  lManager = G4LossTableManager::Instance();
  lManager->Register(this);
  G4LossTableBuilder* bld = lManager->GetTableBuilder();
  theDensityFactor = bld->GetDensityFactors();
  theDensityIdx = bld->GetCoupleIndexes();

  secID = fluoID = augerID = biasID = -1;
  mainSecondaries = 100;
  if("phot" == GetProcessName() || "compt" == GetProcessName()
      || "e-_G4DNAIonisation" == GetProcessName()
      || "hydrogen_G4DNAIonisation" == GetProcessName()
      || "helium_G4DNAIonisation" == GetProcessName()
      || "alpha_G4DNAIonisation" == GetProcessName()
      || "alpha+_G4DNAIonisation" == GetProcessName()
      || "proton_G4DNAIonisation" == GetProcessName()
      || "GenericIon_G4DNAIonisation" == GetProcessName() ) 
    { 
      mainSecondaries = 1; 
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmProcess::~G4VEmProcess()
{
  /*
  if(1 < verboseLevel) {
    G4cout << "G4VEmProcess destruct " << GetProcessName() 
           << "  " << this << "  " <<  theLambdaTable <<G4endl;
  }
  */
  if(isTheMaster) {
    delete theData;
    theData = nullptr;
  }
  delete modelManager;
  delete biasManager;
  lManager->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::Clear()
{
  currentCouple = nullptr;
  preStepLambda = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEmProcess::MinPrimaryEnergy(const G4ParticleDefinition*,
                                        const G4Material*)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::AddEmModel(G4int order, G4VEmModel* p, 
                              const G4Region* region)
{
  G4VEmFluctuationModel* fm = nullptr;
  modelManager->AddEmModel(order, p, fm, region);
  if(p) { p->SetParticleChange(pParticleChange); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::SetEmModel(G4VEmModel* ptr, G4int) 
{
  for(auto & em : emModels) { if(em == ptr) { return; } }
  emModels.push_back(ptr);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmModel* G4VEmProcess::EmModel(size_t index) const
{
  return (index < emModels.size()) ? emModels[index] : nullptr; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::UpdateEmModel(const G4String& nam, 
                                 G4double emin, G4double emax)
{
  modelManager->UpdateEmModel(nam, emin, emax);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4VEmProcess::GetNumberOfModels() const
{
  return modelManager->NumberOfModels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4VEmProcess::GetNumberOfRegionModels(size_t couple_index) const
{
  return modelManager->NumberOfRegionModels(couple_index);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmModel* G4VEmProcess::GetRegionModel(G4int idx, size_t couple_index) const
{
  return modelManager->GetRegionModel(idx, couple_index);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmModel* G4VEmProcess::GetModelByIndex(G4int idx, G4bool ver) const
{
  return modelManager->GetModel(idx, ver);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  isTheMaster = lManager->IsMaster();

  if(!particle) { SetParticle(&part); }

  if(part.GetParticleType() == "nucleus" && 
     part.GetParticleSubType() == "generic") {

    G4String pname = part.GetParticleName();
    if(pname != "deuteron" && pname != "triton" &&
       pname != "alpha" && pname != "He3" &&
       pname != "alpha+"   && pname != "helium" &&
       pname != "hydrogen") {

      particle = G4GenericIon::GenericIon();
      isIon = true;
    }
  }

  if(1 < verboseLevel) {
    G4cout << "G4VEmProcess::PreparePhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << " local particle " << particle->GetParticleName() 
           << G4endl;
  }

  if(particle != &part) { return; }

  G4LossTableBuilder* bld = lManager->GetTableBuilder();

  lManager->PreparePhysicsTable(&part, this, isTheMaster);

  Clear();
  InitialiseProcess(particle);

  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  size_t n = theCoupleTable->GetTableSize();

  theEnergyOfCrossSectionMax.resize(n, 0.0);
  theCrossSectionMax.resize(n, DBL_MAX);

  // initialisation of the process  
  if(!actMinKinEnergy) { minKinEnergy = theParameters->MinKinEnergy(); }
  if(!actMaxKinEnergy) { maxKinEnergy = theParameters->MaxKinEnergy(); }
  if(!actSpline) { splineFlag = theParameters->Spline(); }

  if(isTheMaster) { 
    SetVerboseLevel(theParameters->Verbose());
    if(!theData) { theData = new G4EmDataHandler(2); }
  } else {  
    SetVerboseLevel(theParameters->WorkerVerbose()); 
  }
  applyCuts       = theParameters->ApplyCuts();
  lambdaFactor    = theParameters->LambdaFactor();
  logLambdaFactor = G4Log(lambdaFactor);
  theParameters->DefineRegParamForEM(this);

  // initialisation of models
  numberOfModels = modelManager->NumberOfModels();
  for(G4int i=0; i<numberOfModels; ++i) {
    G4VEmModel* mod = modelManager->GetModel(i);
    if(0 == i) { currentModel = mod; }
    mod->SetPolarAngleLimit(theParameters->MscThetaLimit());
    mod->SetMasterThread(isTheMaster);
    if(mod->HighEnergyLimit() > maxKinEnergy) {
      mod->SetHighEnergyLimit(maxKinEnergy);
    }
  }

  if(lManager->AtomDeexcitation()) { modelManager->SetFluoFlag(true); }
  theCuts = modelManager->Initialise(particle,secondaryParticle,
                                     2.,verboseLevel);
  theCutsGamma    = theCoupleTable->GetEnergyCutsVector(idxG4GammaCut);
  theCutsElectron = theCoupleTable->GetEnergyCutsVector(idxG4ElectronCut);
  theCutsPositron = theCoupleTable->GetEnergyCutsVector(idxG4PositronCut);

  // prepare tables
  if(buildLambdaTable && isTheMaster){
    theLambdaTable = theData->MakeTable(0);
    bld->InitialiseBaseMaterials(theLambdaTable);
  }
  // high energy table
  if(isTheMaster && minKinEnergyPrim < maxKinEnergy){
    theLambdaTablePrim = theData->MakeTable(1);
    bld->InitialiseBaseMaterials(theLambdaTablePrim);
  }
  bld->InitialiseBaseMaterials();
  // forced biasing
  if(biasManager) { 
    biasManager->Initialise(part,GetProcessName(),verboseLevel); 
    biasFlag = false; 
  }
  // defined ID of secondary particles
  G4String nam1 = GetProcessName();
  secID   = G4PhysicsModelCatalog::Register(nam1); 
  if(100 > mainSecondaries) {
    G4String nam2 = nam1 + "_fluo" ;
    G4String nam3 = nam1 + "_auger";
    G4String nam4 = nam1 + "_split";
    fluoID  = G4PhysicsModelCatalog::Register(nam2); 
    augerID = G4PhysicsModelCatalog::Register(nam3); 
    biasID  = G4PhysicsModelCatalog::Register(nam4);
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  if(!masterProc) {
    if(isTheMaster) { masterProc = this; }
    else { masterProc = static_cast<const G4VEmProcess*>(GetMasterProcess());}
  }

  G4String num = part.GetParticleName();
  if(1 < verboseLevel) {
    G4cout << "### G4VEmProcess::BuildPhysicsTable() for "
           << GetProcessName()
           << " and particle " << num
           << " buildLambdaTable= " << buildLambdaTable
           << " isTheMaster= " << isTheMaster 
           << "  " << masterProc 
           << G4endl;
  }

  if(particle == &part) { 

    // worker initialisation
    if(!isTheMaster) {
      theLambdaTable = masterProc->LambdaTable();
      theLambdaTablePrim = masterProc->LambdaTablePrim();

      if(theLambdaTable) { FindLambdaMax(); }

      // local initialisation of models
      G4bool printing = true;
      numberOfModels = modelManager->NumberOfModels();
      for(G4int i=0; i<numberOfModels; ++i) {
        G4VEmModel* mod = GetModelByIndex(i, printing);
        G4VEmModel* mod0= masterProc->GetModelByIndex(i, printing);
        //G4cout << i << ".  " << mod << "   " << mod0 << "  " 
        //     << particle->GetParticleName() << G4endl;
        mod->InitialiseLocal(particle, mod0);
      }
    // master thread
    } else {
      if(buildLambdaTable || minKinEnergyPrim < maxKinEnergy) {
        BuildLambdaTable();
      }
    }
  }

  // explicitly defined printout by particle name
  if(1 < verboseLevel || 
     (0 < verboseLevel && (num == "gamma" || num == "e-" || 
                           num == "e+"    || num == "mu+" || 
                           num == "mu-"   || num == "proton"|| 
                           num == "pi+"   || num == "pi-" || 
                           num == "kaon+" || num == "kaon-" || 
                           num == "alpha" || num == "anti_proton" || 
                           num == "GenericIon"|| num == "alpha++" ||
                           num == "alpha+" || num == "helium" ||
                           num == "hydrogen")))
    { 
      StreamInfo(G4cout, part);
    }

  if(1 < verboseLevel) {
    G4cout << "### G4VEmProcess::BuildPhysicsTable() done for "
           << GetProcessName()
           << " and particle " << num
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::BuildLambdaTable()
{
  if(1 < verboseLevel) {
    G4cout << "G4EmProcess::BuildLambdaTable() for process "
           << GetProcessName() << " and particle "
           << particle->GetParticleName() << "  " << this
           << G4endl;
  }

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  G4LossTableBuilder* bld = lManager->GetTableBuilder();

  G4PhysicsLogVector* aVector = nullptr;
  G4PhysicsLogVector* aVectorPrim = nullptr;
  G4PhysicsLogVector* bVectorPrim = nullptr;

  G4double scale = theParameters->MaxKinEnergy()/theParameters->MinKinEnergy();
  G4int nbin = 
    theParameters->NumberOfBinsPerDecade()*G4lrint(std::log10(scale));
  scale = G4Log(scale);
  if(actBinning) { nbin = std::max(nbin, nLambdaBins); }
  G4double emax1 = std::min(maxKinEnergy, minKinEnergyPrim);
    
  for(size_t i=0; i<numOfCouples; ++i) {

    if (bld->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = 
        theCoupleTable->GetMaterialCutsCouple(i);

      // build main table
      if(buildLambdaTable) {
        delete (*theLambdaTable)[i];

        // if start from zero then change the scale
        G4double emin = minKinEnergy;
        G4bool startNull = false;
        if(startFromNull) {
          G4double e = MinPrimaryEnergy(particle,couple->GetMaterial());
          if(e >= emin) {
            emin = e;
            startNull = true;
          }
        }
        G4double emax = emax1;
        if(emax <= emin) { emax = 2*emin; }
        G4int bin = G4lrint(nbin*G4Log(emax/emin)/scale);
        if(bin < 3) { bin = 3; }
        aVector = new G4PhysicsLogVector(emin, emax, bin);
        aVector->SetSpline(splineFlag);
        modelManager->FillLambdaVector(aVector, couple, startNull);
        if(splineFlag) { aVector->FillSecondDerivatives(); }
        G4PhysicsTableHelper::SetPhysicsVector(theLambdaTable, i, aVector);
      }
      // build high energy table 
      if(minKinEnergyPrim < maxKinEnergy) { 
        delete (*theLambdaTablePrim)[i];

        // start not from zero
        if(!bVectorPrim) {
          G4int bin = G4lrint(nbin*G4Log(maxKinEnergy/minKinEnergyPrim)/scale);
          if(bin < 3) { bin = 3; }
          aVectorPrim = 
            new G4PhysicsLogVector(minKinEnergyPrim, maxKinEnergy, bin);
          bVectorPrim = aVectorPrim;
        } else {
          aVectorPrim = new G4PhysicsLogVector(*bVectorPrim);
        }
        // always use spline
        aVectorPrim->SetSpline(splineFlag);
        modelManager->FillLambdaVector(aVectorPrim, couple, false, 
                                       fIsCrossSectionPrim);
        aVectorPrim->FillSecondDerivatives();
        G4PhysicsTableHelper::SetPhysicsVector(theLambdaTablePrim, i, 
                                               aVectorPrim);
      }
    }
  }

  if(buildLambdaTable) { FindLambdaMax(); }

  if(1 < verboseLevel) {
    G4cout << "Lambda table is built for "
           << particle->GetParticleName()
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::StreamInfo(std::ostream& out, 
                  const G4ParticleDefinition& part, G4bool rst) const
{
  G4String indent = (rst ? "  " : "");
  out << std::setprecision(6);
  out << G4endl << indent << GetProcessName() << ": ";
  if (!rst) {
    out << " for " << part.GetParticleName();
    if (integral) { out << ","; }
  }
  if(integral)  { out << " integral:1 "; }
  if(applyCuts) { out << " applyCuts:1 "; }
  out << " SubType=" << GetProcessSubType();
  if(biasFactor != 1.0) { out << "  BiasingFactor= " << biasFactor; }
  out << " BuildTable=" << buildLambdaTable << G4endl;
  if(buildLambdaTable) {
    if(particle == &part) { 
      size_t length = theLambdaTable->length();
      for(size_t i=0; i<length; ++i) {
        G4PhysicsVector* v = (*theLambdaTable)[i];
        if(v) {
          out << "      Lambda table from ";
          G4double emin = v->Energy(0);
          G4double emax = v->GetMaxEnergy();
          G4int nbin = v->GetVectorLength() - 1;
          if(emin > minKinEnergy) { out << "threshold "; }
          else { out << G4BestUnit(emin,"Energy"); } 
          out << " to "
              << G4BestUnit(emax,"Energy")
              << ", " << G4lrint(nbin/std::log10(emax/emin))
              << " bins/decade, spline: " 
              << splineFlag << G4endl;
          break;
        }
      }
    } else {
      out << "      Used Lambda table of " 
      << particle->GetParticleName() << G4endl;
    }
  }
  if(minKinEnergyPrim < maxKinEnergy) {
    if(particle == &part) { 
      size_t length = theLambdaTablePrim->length();
      for(size_t i=0; i<length; ++i) {
        G4PhysicsVector* v = (*theLambdaTablePrim)[i];
        if(v) { 
          out << "      LambdaPrime table from "
              << G4BestUnit(v->Energy(0),"Energy") 
              << " to "
              << G4BestUnit(v->GetMaxEnergy(),"Energy")
              << " in " << v->GetVectorLength()-1
              << " bins " << G4endl;
          break;
        }
      }
    } else {
      out << "      Used LambdaPrime table of " 
               << particle->GetParticleName() << G4endl;
    }
  }
  StreamProcessInfo(out);
  modelManager->DumpModelList(out, verboseLevel);

  if(verboseLevel > 2 && buildLambdaTable) {
    out << "      LambdaTable address= " << theLambdaTable << G4endl;
    if(theLambdaTable && particle == &part) { 
      out << (*theLambdaTable) << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::StartTracking(G4Track* track)
{
  // reset parameters for the new track
  currentParticle = track->GetParticleDefinition();
  theNumberOfInteractionLengthLeft = -1.0;
  mfpKinEnergy = DBL_MAX;

  if(isIon) { massRatio = proton_mass_c2/currentParticle->GetPDGMass(); }

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

G4double G4VEmProcess::PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition)
{
  *condition = NotForced;
  G4double x = DBL_MAX;

  DefineMaterial(track.GetMaterialCutsCouple());
  preStepKinEnergy      = track.GetKineticEnergy();
  preStepLogKinEnergy   = track.GetDynamicParticle()->GetLogKineticEnergy();
  G4double scaledEnergy = preStepKinEnergy*massRatio;
  SelectModel(scaledEnergy, currentCoupleIndex);

  if(!currentModel->IsActive(scaledEnergy)) { 
    theNumberOfInteractionLengthLeft = -1.0;
    currentInteractionLength = DBL_MAX;
    return x; 
  }
 
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
  if(preStepKinEnergy < mfpKinEnergy) {
    if (integral) {
      ComputeIntegralLambda(preStepKinEnergy, preStepLogKinEnergy);
    } else {
      preStepLambda = GetCurrentLambda(preStepKinEnergy, preStepLogKinEnergy);
    }

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

      theNumberOfInteractionLengthLeft -= 
        previousStepSize/currentInteractionLength;
      theNumberOfInteractionLengthLeft = 
        std::max(theNumberOfInteractionLengthLeft, 0.0);
    }

    // new mean free path and step limit for the next step
    currentInteractionLength = 1.0/preStepLambda;
    x = theNumberOfInteractionLengthLeft * currentInteractionLength;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::ComputeIntegralLambda(G4double e, G4double loge)
{
  // condition to skip recomputation of cross section
  const G4double epeak = theEnergyOfCrossSectionMax[currentCoupleIndex];
  if(e <= epeak && e/lambdaFactor >= mfpKinEnergy) { return; }

  // recomputation is needed 
  if (e <= epeak) {
    preStepLambda = GetCurrentLambda(e, loge);
    mfpKinEnergy  = e;
  } else {
    const G4double e1 = e*lambdaFactor;
    if (e1 > epeak) {
      preStepLambda = GetCurrentLambda(e, loge);
      mfpKinEnergy  = e;
      const G4double preStepLambda1 = GetCurrentLambda(e1,loge+logLambdaFactor);
      if (preStepLambda1 > preStepLambda) {
        mfpKinEnergy  = e1;
        preStepLambda = preStepLambda1;
      }
    } else {
      preStepLambda = fFactor*theCrossSectionMax[currentCoupleIndex];
      mfpKinEnergy  = epeak;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEmProcess::PostStepDoIt(const G4Track& track,
                                              const G4Step& step)
{
  // In all cases clear number of interaction lengths
  theNumberOfInteractionLengthLeft = -1.0;
  mfpKinEnergy = DBL_MAX; 

  fParticleChange.InitializeForPostStep(track);

  // Do not make anything if particle is stopped, the annihilation then
  // should be performed by the AtRestDoIt!
  if (track.GetTrackStatus() == fStopButAlive) { return &fParticleChange; }

  const G4double finalT    = track.GetKineticEnergy();
  const G4double logFinalT = track.GetDynamicParticle()->GetLogKineticEnergy();

  // forced process - should happen only once per track
  if(biasFlag) {
    if(biasManager->ForcedInteractionRegion(currentCoupleIndex)) {
      biasFlag = false;
    }
  }

  // Integral approach
  if (integral) {
    G4double lx = GetLambda(finalT, currentCouple, logFinalT);
    if(preStepLambda<lx && 1 < verboseLevel) {
      G4cout << "WARNING: for " << currentParticle->GetParticleName() 
             << " and " << GetProcessName()
             << " E(MeV)= " << finalT/MeV
             << " preLambda= " << preStepLambda << " < " 
             << lx << " (postLambda) "
             << G4endl;  
    }

    if(preStepLambda*G4UniformRand() > lx) {
      ClearNumberOfInteractionLengthLeft();
      return &fParticleChange;
    }
  }

  G4double scaledEnergy = finalT*massRatio;
  SelectModel(scaledEnergy, currentCoupleIndex);
  if(!currentModel->IsActive(scaledEnergy)) { return &fParticleChange; }

  // define new weight for primary and secondaries
  G4double weight = fParticleChange.GetParentWeight();
  if(weightFlag) { 
    weight /= biasFactor; 
    fParticleChange.ProposeWeight(weight);
  }
  

  if(1 < verboseLevel) {
    G4cout << "G4VEmProcess::PostStepDoIt: Sample secondary; E= "
           << finalT/MeV
           << " MeV; model= (" << currentModel->LowEnergyLimit()
           << ", " <<  currentModel->HighEnergyLimit() << ")"
           << G4endl;
  }


  // sample secondaries
  secParticles.clear();
  currentModel->SampleSecondaries(&secParticles, 
                                  currentCouple, 
                                  track.GetDynamicParticle(),
                                  (*theCuts)[currentCoupleIndex]);

  G4int num0 = secParticles.size();

  // splitting or Russian roulette
  if(biasManager) {
    if(biasManager->SecondaryBiasingRegion(currentCoupleIndex)) {
      G4double eloss = 0.0;
      weight *= biasManager->ApplySecondaryBiasing(
        secParticles, track, currentModel, &fParticleChange, eloss, 
        currentCoupleIndex, (*theCuts)[currentCoupleIndex],
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
    G4double edep = fParticleChange.GetLocalEnergyDeposit();
    G4double time = track.GetGlobalTime();
     
    for (G4int i=0; i<num; ++i) {
      if (secParticles[i]) {
        G4DynamicParticle* dp = secParticles[i];
        const G4ParticleDefinition* p = dp->GetParticleDefinition();
        G4double e = dp->GetKineticEnergy();
        G4bool good = true;
        if(applyCuts) {
          if (p == theGamma) {
            if (e < (*theCutsGamma)[currentCoupleIndex]) { good = false; }

          } else if (p == theElectron) {
            if (e < (*theCutsElectron)[currentCoupleIndex]) { good = false; }

          } else if (p == thePositron) {
            if (electron_mass_c2 < (*theCutsGamma)[currentCoupleIndex] &&
                e < (*theCutsPositron)[currentCoupleIndex]) {
              good = false;
              e += 2.0*electron_mass_c2;
            }
          }
          // added secondary if it is good
        }
        if (good) { 
          G4Track* t = new G4Track(dp, time, track.GetPosition());
          t->SetTouchableHandle(track.GetTouchableHandle());
          if (biasManager) {
            t->SetWeight(weight * biasManager->GetWeight(i));
          } else {
            t->SetWeight(weight);
          }
          pParticleChange->AddSecondary(t);

          // define type of secondary
          if(i < mainSecondaries) { t->SetCreatorModelIndex(secID); }
          else if(i < num0) {
            if(p == theGamma) { 
              t->SetCreatorModelIndex(fluoID);
            } else {
              t->SetCreatorModelIndex(augerID);
            }
          } else {
            t->SetCreatorModelIndex(biasID);
          }
          /* 
          G4cout << "Secondary(post step) has weight " << t->GetWeight() 
                 << ", Ekin= " << t->GetKineticEnergy()/MeV << " MeV "
                 << GetProcessName() << " fluoID= " << fluoID
                 << " augerID= " << augerID <<G4endl;
          */
        } else {
          delete dp;
          edep += e;
        }
      } 
    }
    fParticleChange.ProposeLocalEnergyDeposit(edep);
  }

  if(0.0 == fParticleChange.GetProposedKineticEnergy() &&
     fAlive == fParticleChange.GetTrackStatus()) {
    if(particle->GetProcessManager()->GetAtRestProcessVector()->size() > 0)
         { fParticleChange.ProposeTrackStatus(fStopButAlive); }
    else { fParticleChange.ProposeTrackStatus(fStopAndKill); }
  }

  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4VEmProcess::StorePhysicsTable(const G4ParticleDefinition* part,
                                       const G4String& directory,
                                       G4bool ascii)
{
  G4bool yes = true;
  if(!isTheMaster) { return yes; }

  if ( theLambdaTable && part == particle) {
    const G4String& nam = 
      GetPhysicsTableFileName(part,directory,"Lambda",ascii);
    yes = theLambdaTable->StorePhysicsTable(nam,ascii);

    if ( yes ) {
      if(0 < verboseLevel) G4cout << "Stored: " << nam << G4endl;
    } else {
      G4cout << "Fail to store Physics Table for " 
             << particle->GetParticleName()
             << " and process " << GetProcessName()
             << " in the directory <" << directory
             << "> " << G4endl;
    }
  }
  if ( theLambdaTablePrim && part == particle) {
    const G4String& name = 
      GetPhysicsTableFileName(part,directory,"LambdaPrim",ascii);
    yes = theLambdaTablePrim->StorePhysicsTable(name,ascii);

    if ( yes ) {
      if(0 < verboseLevel) {
	G4cout << "Physics table prim is stored for " 
	       << particle->GetParticleName()
	       << " and process " << GetProcessName()
	       << " in the directory <" << directory
	       << "> " << G4endl;
      }
    } else {
      G4cout << "Fail to store Physics Table Prim for " 
             << particle->GetParticleName()
             << " and process " << GetProcessName()
             << " in the directory <" << directory
             << "> " << G4endl;
    }
  }
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool G4VEmProcess::RetrievePhysicsTable(const G4ParticleDefinition* part,
                                          const G4String& directory,
                                          G4bool ascii)
{
  if(1 < verboseLevel) {
    G4cout << "G4VEmProcess::RetrievePhysicsTable() for "
           << part->GetParticleName() << " and process "
           << GetProcessName() << G4endl;
  }
  G4bool yes = true;

  if((!buildLambdaTable && minKinEnergyPrim > maxKinEnergy) 
     || particle != part) { return yes; }

  const G4String particleName = part->GetParticleName();

  if(buildLambdaTable) {
    const G4String& filename = 
      GetPhysicsTableFileName(part,directory,"Lambda",ascii);
    yes = G4PhysicsTableHelper::RetrievePhysicsTable(theLambdaTable,
                                                     filename,ascii);
    if ( yes ) {
      if (0 < verboseLevel) {
        G4cout << "Lambda table for " << particleName 
               << " is Retrieved from <"
               << filename << ">"
               << G4endl;
      }
      if(theParameters->Spline()) {
        size_t n = theLambdaTable->length();
        for(size_t i=0; i<n; ++i) {
          if((* theLambdaTable)[i]) {
            (* theLambdaTable)[i]->SetSpline(true);
          }
        }
      }
    } else {
      if (1 < verboseLevel) {
        G4cout << "Lambda table for " << particleName << " in file <"
               << filename << "> is not exist"
               << G4endl;
      }
    }
  }
  if(minKinEnergyPrim < maxKinEnergy) {
    const G4String& filename = 
      GetPhysicsTableFileName(part,directory,"LambdaPrim",ascii);
    yes = G4PhysicsTableHelper::RetrievePhysicsTable(theLambdaTablePrim,
                                                     filename,ascii);
    if ( yes ) {
      if (0 < verboseLevel) {
        G4cout << "Lambda table prim for " << particleName 
               << " is Retrieved from <"
               << filename << ">"
               << G4endl;
      }
      if(theParameters->Spline()) {
        size_t n = theLambdaTablePrim->length();
        for(size_t i=0; i<n; ++i) {
          if((* theLambdaTablePrim)[i]) {
            (* theLambdaTablePrim)[i]->SetSpline(true);
          }
        }
      }
    } else {
      if (1 < verboseLevel) {
        G4cout << "Lambda table prim for " << particleName << " in file <"
               << filename << "> is not exist"
               << G4endl;
      }
    }
  }
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4VEmProcess::CrossSectionPerVolume(G4double kineticEnergy,
                                    const G4MaterialCutsCouple* couple,
                                    G4double logKinEnergy)
{
  // Cross section per atom is calculated
  DefineMaterial(couple);
  G4double cross = 0.0;
  if(buildLambdaTable) {
    cross = GetCurrentLambda(kineticEnergy, 
      (logKinEnergy < DBL_MAX) ? logKinEnergy : G4Log(kineticEnergy));
  } else {
    SelectModel(kineticEnergy, currentCoupleIndex);
    if(currentModel) {
      cross = fFactor*currentModel->CrossSectionPerVolume(currentMaterial,
                                                          currentParticle,
                                                          kineticEnergy);
    }
  }
  return std::max(cross, 0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEmProcess::GetMeanFreePath(const G4Track& track,
                                       G4double,
                                       G4ForceCondition* condition)
{
  *condition = NotForced;
  return G4VEmProcess::MeanFreePath(track);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEmProcess::MeanFreePath(const G4Track& track)
{
  const G4double kinEnergy = track.GetKineticEnergy();
  CurrentSetup(track.GetMaterialCutsCouple(), kinEnergy);
  const G4double xs = GetCurrentLambda(kinEnergy,
                             track.GetDynamicParticle()->GetLogKineticEnergy());
  return (0.0 < xs) ? 1.0/xs : DBL_MAX; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4VEmProcess::ComputeCrossSectionPerAtom(G4double kinEnergy, 
                                         G4double Z, G4double A, G4double cut)
{
  SelectModel(kinEnergy, currentCoupleIndex);
  return (currentModel) ? 
    currentModel->ComputeCrossSectionPerAtom(currentParticle, kinEnergy,
                                             Z, A, cut) : 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::FindLambdaMax()
{
  if(1 < verboseLevel) {
    G4cout << "### G4VEmProcess::FindLambdaMax: " 
           << particle->GetParticleName() 
           << " and process " << GetProcessName() << "  " << G4endl; 
  }
  size_t n = theLambdaTable->length();
  G4PhysicsVector* pv;
  G4double e, ss, emax, smax;

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* 
G4VEmProcess::LambdaPhysicsVector(const G4MaterialCutsCouple*)
{
  G4PhysicsVector* v = 
    new G4PhysicsLogVector(minKinEnergy, maxKinEnergy, nLambdaBins);
  v->SetSpline(theParameters->Spline());
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Element* G4VEmProcess::GetCurrentElement() const
{
  const G4Element* elm = 
    (currentModel) ? currentModel->GetCurrentElement() : nullptr; 
  return elm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::SetCrossSectionBiasingFactor(G4double f, G4bool flag)
{
  if(f > 0.0) { 
    biasFactor = f; 
    weightFlag = flag;
    if(1 < verboseLevel) {
      G4cout << "### SetCrossSectionBiasingFactor: for " 
             << particle->GetParticleName() 
             << " and process " << GetProcessName()
             << " biasFactor= " << f << " weightFlag= " << flag 
             << G4endl; 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4VEmProcess::ActivateForcedInteraction(G4double length, const G4String& r,
                                        G4bool flag)
{
  if(!biasManager) { biasManager = new G4EmBiasingManager(); }
  if(1 < verboseLevel) {
    G4cout << "### ActivateForcedInteraction: for " 
           << particle->GetParticleName() 
           << " and process " << GetProcessName()
           << " length(mm)= " << length/mm
           << " in G4Region <" << r 
           << "> weightFlag= " << flag 
           << G4endl; 
  }
  weightFlag = flag;
  biasManager->ActivateForcedInteraction(length, r);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4VEmProcess::ActivateSecondaryBiasing(const G4String& region,
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

void G4VEmProcess::SetLambdaBinning(G4int n)
{
  if(5 < n && n < 10000000) {  
    nLambdaBins = n; 
    actBinning = true;
  } else { 
    G4double e = (G4double)n;
    PrintWarning("SetLambdaBinning", e); 
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::SetMinKinEnergy(G4double e)
{
  if(1.e-3*eV < e && e < maxKinEnergy) { 
    nLambdaBins = G4lrint(nLambdaBins*G4Log(maxKinEnergy/e)
                          /G4Log(maxKinEnergy/minKinEnergy));
    minKinEnergy = e;
    actMinKinEnergy = true;
  } else { PrintWarning("SetMinKinEnergy", e); } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::SetMaxKinEnergy(G4double e)
{
  if(minKinEnergy < e && e < 1.e+6*TeV) { 
    nLambdaBins = G4lrint(nLambdaBins*G4Log(e/minKinEnergy)
                          /G4Log(maxKinEnergy/minKinEnergy));
    maxKinEnergy = e;
    actMaxKinEnergy = true;
  } else { PrintWarning("SetMaxKinEnergy", e); } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::SetMinKinEnergyPrim(G4double e)
{
  if(theParameters->MinKinEnergy() <= e && 
     e <= theParameters->MaxKinEnergy()) { minKinEnergyPrim = e; } 
  else { PrintWarning("SetMinKinEnergyPrim", e); } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmProcess* G4VEmProcess::GetEmProcess(const G4String& nam)
{
  return (nam == GetProcessName()) ? this : nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4VEmProcess::GetLambda(G4double kinEnergy, const G4MaterialCutsCouple* couple)
{
  CurrentSetup(couple, kinEnergy);
  return GetCurrentLambda(kinEnergy, G4Log(kinEnergy));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::PrintWarning(G4String tit, G4double val)
{
  G4String ss = "G4VEmProcess::" + tit;
  G4ExceptionDescription ed;
  ed << "Parameter is out of range: " << val 
     << " it will have no effect!\n" << "  Process " 
     << GetProcessName() << "  nbins= " << theParameters->NumberOfBins()
     << " Emin(keV)= " << theParameters->MinKinEnergy()/keV 
     << " Emax(GeV)= " << theParameters->MaxKinEnergy()/GeV;
  G4Exception(ss, "em0044", JustWarning, ed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::ProcessDescription(std::ostream& out) const
{
  if(particle) {
    StreamInfo(out, *particle, true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
