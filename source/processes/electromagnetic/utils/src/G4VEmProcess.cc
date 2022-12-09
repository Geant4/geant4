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
#include "G4EmParameters.hh"
#include "G4EmProcessSubType.hh"
#include "G4EmTableUtil.hh"
#include "G4EmUtility.hh"
#include "G4DNAModelSubType.hh"
#include "G4GenericIon.hh"
#include "G4Log.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmProcess::G4VEmProcess(const G4String& name, G4ProcessType type):
  G4VDiscreteProcess(name, type)
{
  theParameters = G4EmParameters::Instance();
  SetVerboseLevel(1);

  // Size of tables
  minKinEnergy = 0.1*CLHEP::keV;
  maxKinEnergy = 100.0*CLHEP::TeV;

  // default lambda factor
  invLambdaFactor = 1.0/lambdaFactor;

  // particle types
  theGamma = G4Gamma::Gamma();
  theElectron = G4Electron::Electron();
  thePositron = G4Positron::Positron();

  pParticleChange = &fParticleChange;
  fParticleChange.SetSecondaryWeightByProcess(true);
  secParticles.reserve(5);

  modelManager = new G4EmModelManager();
  lManager = G4LossTableManager::Instance();
  lManager->Register(this);
  isTheMaster = lManager->IsMaster();
  G4LossTableBuilder* bld = lManager->GetTableBuilder();
  theDensityFactor = bld->GetDensityFactors();
  theDensityIdx = bld->GetCoupleIndexes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmProcess::~G4VEmProcess()
{
  if(isTheMaster) {
    delete theData;
    delete theEnergyOfCrossSectionMax;
  }
  delete modelManager;
  delete biasManager;
  lManager->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::AddEmModel(G4int order, G4VEmModel* ptr, 
                              const G4Region* region)
{
  if(nullptr == ptr) { return; }
  G4VEmFluctuationModel* fm = nullptr;
  modelManager->AddEmModel(order, ptr, fm, region);
  ptr->SetParticleChange(pParticleChange);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::SetEmModel(G4VEmModel* ptr, G4int) 
{
  if(nullptr == ptr) { return; }
  if(!emModels.empty()) {
    for(auto & em : emModels) { if(em == ptr) { return; } }
  }
  emModels.push_back(ptr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  if(nullptr == particle) { SetParticle(&part); }

  if(part.GetParticleType() == "nucleus" && 
     part.GetParticleSubType() == "generic") {

    G4String pname = part.GetParticleName();
    if(pname != "deuteron" && pname != "triton" &&
       pname != "alpha" && pname != "alpha+" &&
       pname != "helium" && pname != "hydrogen") {

      particle = G4GenericIon::GenericIon();
      isIon = true;
    }
  }
  if(particle != &part) { return; }

  lManager->PreparePhysicsTable(&part, this, isTheMaster);

  // for new run
  currentCouple = nullptr;
  preStepLambda = 0.0;
  fLambdaEnergy = 0.0;

  InitialiseProcess(particle);

  G4LossTableBuilder* bld = lManager->GetTableBuilder();
  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  theCutsGamma    = theCoupleTable->GetEnergyCutsVector(idxG4GammaCut);
  theCutsElectron = theCoupleTable->GetEnergyCutsVector(idxG4ElectronCut);
  theCutsPositron = theCoupleTable->GetEnergyCutsVector(idxG4PositronCut);

  // initialisation of the process  
  if(!actMinKinEnergy) { minKinEnergy = theParameters->MinKinEnergy(); }
  if(!actMaxKinEnergy) { maxKinEnergy = theParameters->MaxKinEnergy(); }

  applyCuts       = theParameters->ApplyCuts();
  lambdaFactor    = theParameters->LambdaFactor();
  invLambdaFactor = 1.0/lambdaFactor;
  theParameters->DefineRegParamForEM(this);

  // integral option may be disabled
  if(!theParameters->Integral()) { fXSType = fEmNoIntegral; }

  // prepare tables
  if(isTheMaster) {
    if(nullptr == theData) { theData = new G4EmDataHandler(2); }

    if(buildLambdaTable) {
      theLambdaTable = theData->MakeTable(0);
      bld->InitialiseBaseMaterials(theLambdaTable);
    }
    // high energy table
    if(minKinEnergyPrim < maxKinEnergy) {
      theLambdaTablePrim = theData->MakeTable(1);
      bld->InitialiseBaseMaterials(theLambdaTablePrim);
    }
  }
  // models
  baseMat = bld->GetBaseMaterialFlag();
  numberOfModels = modelManager->NumberOfModels();
  currentModel = modelManager->GetModel(0);
  if(nullptr != lManager->AtomDeexcitation()) { 
    modelManager->SetFluoFlag(true); 
  }
  // forced biasing
  if(nullptr != biasManager) { 
    biasManager->Initialise(part, GetProcessName(), verboseLevel); 
    biasFlag = false;
  }

  theCuts =
    G4EmTableUtil::PrepareEmProcess(this, particle, secondaryParticle,
                                    modelManager, maxKinEnergy, 
                                    secID, tripletID, mainSecondaries,
                                    verboseLevel, isTheMaster);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  if(nullptr == masterProc) {
    if(isTheMaster) { masterProc = this; }
    else { masterProc = static_cast<const G4VEmProcess*>(GetMasterProcess());}
  }
  G4int nModels = modelManager->NumberOfModels();
  G4bool isLocked = theParameters->IsPrintLocked();
  G4bool toBuild = (buildLambdaTable || minKinEnergyPrim < maxKinEnergy);

  G4EmTableUtil::BuildEmProcess(this, masterProc, particle, &part,
                                nModels, verboseLevel, isTheMaster,
                                isLocked, toBuild, baseMat);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::BuildLambdaTable()
{
  G4double scale = theParameters->MaxKinEnergy()/theParameters->MinKinEnergy();
  G4int nbin = 
    theParameters->NumberOfBinsPerDecade()*G4lrint(std::log10(scale));
  if(actBinning) { nbin = std::max(nbin, nLambdaBins); }
  scale = nbin/G4Log(scale);
  
  G4LossTableBuilder* bld = lManager->GetTableBuilder();
  G4EmTableUtil::BuildLambdaTable(this, particle, modelManager,
                                  bld, theLambdaTable, theLambdaTablePrim,
                                  minKinEnergy, minKinEnergyPrim,
                                  maxKinEnergy, scale, verboseLevel,
                                  startFromNull, splineFlag);
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
  }
  if(fXSType != fEmNoIntegral)  { out << " XStype:" << fXSType; }
  if(applyCuts) { out << " applyCuts:1 "; }
  out << " SubType=" << GetProcessSubType();
  if(biasFactor != 1.0) { out << "  BiasingFactor= " << biasFactor; }
  out << " BuildTable=" << buildLambdaTable << G4endl;
  if(buildLambdaTable) {
    if(particle == &part) { 
      for(auto & v : *theLambdaTable) {
        if(nullptr != v) {
          out << "      Lambda table from ";
          G4double emin = v->Energy(0);
          G4double emax = v->GetMaxEnergy();
          G4int nbin = G4int(v->GetVectorLength() - 1);
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
      for(auto & v : *theLambdaTablePrim) {
        if(nullptr != v) {
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
  preStepKinEnergy    = track.GetKineticEnergy();
  preStepLogKinEnergy = track.GetDynamicParticle()->GetLogKineticEnergy();
  const G4double scaledEnergy = preStepKinEnergy*massRatio;
  SelectModel(scaledEnergy, currentCoupleIndex);
  /*
  G4cout << "PostStepGetPhysicalInteractionLength: idx= " << currentCoupleIndex
         << "  couple: " << currentCouple << G4endl;
  */
  if(!currentModel->IsActive(scaledEnergy)) { 
    theNumberOfInteractionLengthLeft = -1.0;
    currentInteractionLength = DBL_MAX;
    return x; 
  }
 
  // forced biasing only for primary particles
  if(biasManager) {
    if(0 == track.GetParentID()) {
      if(biasFlag && 
         biasManager->ForcedInteractionRegion((G4int)currentCoupleIndex)) {
        return biasManager->GetStepLimit((G4int)currentCoupleIndex, previousStepSize);
      }
    }
  }

  // compute mean free path
  ComputeIntegralLambda(preStepKinEnergy, preStepLogKinEnergy);

  // zero cross section
  if(preStepLambda <= 0.0) { 
    theNumberOfInteractionLengthLeft = -1.0;
    currentInteractionLength = DBL_MAX;

  } else {

    // non-zero cross section
    if (theNumberOfInteractionLengthLeft < 0.0) {

      // beggining of tracking (or just after DoIt of this process)
      theNumberOfInteractionLengthLeft = -G4Log( G4UniformRand() );
      theInitialNumberOfInteractionLength = theNumberOfInteractionLengthLeft; 

    } else {

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
  if(fXSType == fEmNoIntegral) {
    preStepLambda = GetCurrentLambda(e, loge);

  } else if(fXSType == fEmIncreasing) {
    if(e*invLambdaFactor < mfpKinEnergy) {
      mfpKinEnergy = e;
      preStepLambda = GetCurrentLambda(e, loge); 
    }

  } else if(fXSType == fEmDecreasing) {
    if(e < mfpKinEnergy) { 
      const G4double e1 = e*lambdaFactor;
      preStepLambda = GetCurrentLambda(e1); 
      mfpKinEnergy = e1;
    }

  } else if(fXSType == fEmOnePeak) {
    const G4double epeak = (*theEnergyOfCrossSectionMax)[currentCoupleIndex];
    if(e <= epeak) {
      if(e*invLambdaFactor < mfpKinEnergy) {
        mfpKinEnergy = e;
        preStepLambda = GetCurrentLambda(e, loge); 
      }
    } else if(e < mfpKinEnergy) { 
      const G4double e1 = std::max(epeak, e*lambdaFactor);
      preStepLambda = GetCurrentLambda(e1); 
      mfpKinEnergy = e1;
    }

  } else {
    preStepLambda = GetCurrentLambda(e, loge); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEmProcess::PostStepDoIt(const G4Track& track,
                                              const G4Step& step)
{
  // clear number of interaction lengths in any case
  theNumberOfInteractionLengthLeft = -1.0;
  mfpKinEnergy = DBL_MAX;

  fParticleChange.InitializeForPostStep(track);

  // Do not make anything if particle is stopped, the annihilation then
  // should be performed by the AtRestDoIt!
  if (track.GetTrackStatus() == fStopButAlive) { return &fParticleChange; }

  const G4double finalT = track.GetKineticEnergy();

  // forced process - should happen only once per track
  if(biasFlag) {
    if(biasManager->ForcedInteractionRegion((G4int)currentCoupleIndex)) {
      biasFlag = false;
    }
  }

  // check active and select model
  const G4double scaledEnergy = finalT*massRatio;
  SelectModel(scaledEnergy, currentCoupleIndex);
  if(!currentModel->IsActive(scaledEnergy)) { return &fParticleChange; }

  // Integral approach
  if (fXSType != fEmNoIntegral) {
    const G4double logFinalT = 
      track.GetDynamicParticle()->GetLogKineticEnergy();
    const G4double lx = std::max(GetCurrentLambda(finalT, logFinalT), 0.0);
#ifdef G4VERBOSE
    if(preStepLambda < lx && 1 < verboseLevel) {
      G4cout << "WARNING: for " << currentParticle->GetParticleName() 
             << " and " << GetProcessName() << " E(MeV)= " << finalT/MeV
             << " preLambda= " << preStepLambda 
             << " < " << lx << " (postLambda) " << G4endl;  
    }
#endif
    // if false interaction then use new cross section value
    // if both values are zero - no interaction
    if(preStepLambda*G4UniformRand() >= lx) {
      return &fParticleChange;
    }
  }

  // define new weight for primary and secondaries
  G4double weight = fParticleChange.GetParentWeight();
  if(weightFlag) { 
    weight /= biasFactor; 
    fParticleChange.ProposeWeight(weight);
  }
  
#ifdef G4VERBOSE
  if(1 < verboseLevel) {
    G4cout << "G4VEmProcess::PostStepDoIt: Sample secondary; E= "
           << finalT/MeV
           << " MeV; model= (" << currentModel->LowEnergyLimit()
           << ", " <<  currentModel->HighEnergyLimit() << ")"
           << G4endl;
  }
#endif

  // sample secondaries
  secParticles.clear();
  currentModel->SampleSecondaries(&secParticles, 
                                  currentCouple, 
                                  track.GetDynamicParticle(),
                                  (*theCuts)[currentCoupleIndex]);

  G4int num0 = (G4int)secParticles.size();

  // splitting or Russian roulette
  if(biasManager) {
    if(biasManager->SecondaryBiasingRegion((G4int)currentCoupleIndex)) {
      G4double eloss = 0.0;
      weight *= biasManager->ApplySecondaryBiasing(
        secParticles, track, currentModel, &fParticleChange, eloss, 
        (G4int)currentCoupleIndex, (*theCuts)[currentCoupleIndex],
        step.GetPostStepPoint()->GetSafety());
      if(eloss > 0.0) {
        eloss += fParticleChange.GetLocalEnergyDeposit();
        fParticleChange.ProposeLocalEnergyDeposit(eloss);
      }
    }
  }

  // save secondaries
  G4int num = (G4int)secParticles.size();
  if(num > 0) {

    fParticleChange.SetNumberOfSecondaries(num);
    G4double edep = fParticleChange.GetLocalEnergyDeposit();
    G4double time = track.GetGlobalTime();

    G4int n1(0), n2(0);
    if(num > mainSecondaries) { 
      currentModel->FillNumberOfSecondaries(n1, n2);
    }
     
    for (G4int i=0; i<num; ++i) {
      G4DynamicParticle* dp = secParticles[i];
      if (nullptr != dp) {
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
          if(i < mainSecondaries) { 
            t->SetCreatorModelID(secID);
            if(GetProcessSubType() == fComptonScattering && p == theGamma) {
              t->SetCreatorModelID(_ComptonGamma);
            }
          } else if(i < mainSecondaries + n1) {
            t->SetCreatorModelID(tripletID);
          } else if(i < mainSecondaries + n1 + n2) {
            t->SetCreatorModelID(_IonRecoil);
          } else {
            if(i < num0) {
              if(p == theGamma) { 
                t->SetCreatorModelID(fluoID);
              } else {
                t->SetCreatorModelID(augerID);
              }
            } else {
              t->SetCreatorModelID(secID);
            }
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
  if(!isTheMaster || part != particle) { return true; }
  if(G4EmTableUtil::StoreTable(this, part, theLambdaTable,
			       directory, "Lambda",
                               verboseLevel, ascii) &&
     G4EmTableUtil::StoreTable(this, part, theLambdaTablePrim,
			       directory, "LambdaPrim",
                               verboseLevel, ascii)) { 
     return true;
  }
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool G4VEmProcess::RetrievePhysicsTable(const G4ParticleDefinition* part,
                                          const G4String& dir,
                                          G4bool ascii)
{
  if(!isTheMaster || part != particle) { return true; }
  G4bool yes = true;
  if(buildLambdaTable) {
    yes = G4EmTableUtil::RetrieveTable(this, part, theLambdaTable, dir,
                                       "Lambda", verboseLevel,
                                       ascii, splineFlag);
  }
  if(yes && minKinEnergyPrim < maxKinEnergy) {
    yes = G4EmTableUtil::RetrieveTable(this, part, theLambdaTablePrim, dir,
                                       "LambdaPrim", verboseLevel,
                                       ascii, splineFlag);
  }
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEmProcess::GetCrossSection(G4double kinEnergy,
                                       const G4MaterialCutsCouple* couple)
{
  CurrentSetup(couple, kinEnergy);
  return GetCurrentLambda(kinEnergy, G4Log(kinEnergy));
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

G4PhysicsVector* 
G4VEmProcess::LambdaPhysicsVector(const G4MaterialCutsCouple* couple)
{
  DefineMaterial(couple);
  G4PhysicsVector* newv = new G4PhysicsLogVector(minKinEnergy, maxKinEnergy, 
                                                 nLambdaBins, splineFlag);
  return newv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Element* G4VEmProcess::GetCurrentElement() const
{
  return (nullptr != currentModel) ?
    currentModel->GetCurrentElement(currentMaterial) : nullptr; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Element* G4VEmProcess::GetTargetElement() const
{
  return (nullptr != currentModel) ?
    currentModel->GetCurrentElement(currentMaterial) : nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Isotope* G4VEmProcess::GetTargetIsotope() const
{
  return (nullptr != currentModel) ?
    currentModel->GetCurrentIsotope(GetCurrentElement()) : nullptr;
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
  if(nullptr == biasManager) { biasManager = new G4EmBiasingManager(); }
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

G4double G4VEmProcess::PolarAngleLimit() const
{
  return theParameters->MscThetaLimit();
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
  if(nullptr != particle) {
    StreamInfo(out, *particle, true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
