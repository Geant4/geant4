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
// Geant4 class G4EmTableUtil
//
// Author V.Ivanchenko 14.03.2022
//

#include "G4EmTableUtil.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCutsTable.hh"
#include "G4EmParameters.hh"
#include "G4EmUtility.hh"
#include "G4LossTableManager.hh"
#include "G4EmTableType.hh"
#include "G4PhysicsModelCatalog.hh"
#include "G4LowEnergyEmProcessSubType.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ProcessManager.hh"
#include "G4UIcommand.hh"
#include "G4GenericIon.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const G4DataVector*
G4EmTableUtil::PrepareEmProcess(G4VEmProcess* proc,
                                const G4ParticleDefinition* part,
                                const G4ParticleDefinition* secPart,
			        G4EmModelManager* modelManager,
                                const G4double& maxKinEnergy,
			        G4int& secID, G4int& tripletID,
                                G4int& mainSec, const G4int& verb,
                                const G4bool& master)
{
  G4EmParameters* param = G4EmParameters::Instance();

  // initialisation of models
  G4double plimit = param->MscThetaLimit();
  G4int nModels = modelManager->NumberOfModels();
  for(G4int i=0; i<nModels; ++i) {
    G4VEmModel* mod = modelManager->GetModel(i);
    if(nullptr == mod) { continue; }
    mod->SetPolarAngleLimit(plimit);
    mod->SetMasterThread(master);
    if(mod->HighEnergyLimit() > maxKinEnergy) {
      mod->SetHighEnergyLimit(maxKinEnergy);
    }
    proc->SetEmModel(mod);
  }

  // defined ID of secondary particles and verbosity
  G4int stype = proc->GetProcessSubType();
  if(stype == fAnnihilation) {
    secID = _Annihilation;
    tripletID = _TripletGamma;
  } else if(stype == fGammaConversion) {
    secID = _PairProduction;
    mainSec = 2;
  } else if(stype == fPhotoElectricEffect) {
    secID = _PhotoElectron;
  } else if(stype == fComptonScattering) {
    secID = _ComptonElectron;
  } else if(stype >= fLowEnergyElastic) {
    secID = fDNAUnknownModel;
  }
  if(master) { 
    proc->SetVerboseLevel(param->Verbose());
  } else {  
    proc->SetVerboseLevel(param->WorkerVerbose()); 
  }

  // model initialisation
  const G4DataVector* cuts = modelManager->Initialise(part, secPart, verb);

  if(1 < verb) {
    G4cout << "### G4VEmProcess::PreparePhysicsTable() done for " 
           << proc->GetProcessName()
           << " and particle " << part->GetParticleName()
           << G4endl;
  }
  return cuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4EmTableUtil::BuildEmProcess(G4VEmProcess* proc,
                                   const G4VEmProcess* masterProc,
			           const G4ParticleDefinition* firstPart,
			           const G4ParticleDefinition* part,
			           const G4int nModels, const G4int verb,
                                   const G4bool master, const G4bool isLocked,
                                   const G4bool toBuild, G4bool& baseMat)
{
  G4String num = part->GetParticleName();
  if(1 < verb) {
    G4cout << "### G4VEmProcess::BuildPhysicsTable() for "
           << proc->GetProcessName() << " and particle " << num
           << " buildLambdaTable=" << toBuild << " master= " << master 
           << G4endl;
  }

  if(firstPart == part) { 

    // worker initialisation
    if(!master) { 
      proc->SetLambdaTable(masterProc->LambdaTable());
      proc->SetLambdaTablePrim(masterProc->LambdaTablePrim());
      proc->SetCrossSectionType(masterProc->CrossSectionType());
      proc->SetEnergyOfCrossSectionMax(masterProc->EnergyOfCrossSectionMax());

      // local initialisation of models
      baseMat = masterProc->UseBaseMaterial();
      G4bool printing = true;
      for(G4int i=0; i<nModels; ++i) {
        G4VEmModel* mod = proc->GetModelByIndex(i, printing);
        G4VEmModel* mod0= masterProc->GetModelByIndex(i, printing);
        mod->SetUseBaseMaterials(baseMat);
        mod->InitialiseLocal(part, mod0);
      }
      // master thread
    } else {
      if(toBuild) { proc->BuildLambdaTable(); }
      auto fXSType = proc->CrossSectionType();
      auto v = proc->EnergyOfCrossSectionMax();
      delete v;
      v = nullptr;
      if(fXSType == fEmOnePeak) {
        auto table = proc->LambdaTable();
        if(nullptr == table) {
	  v = G4EmUtility::FindCrossSectionMax(proc, part);
	} else {
	  v = G4EmUtility::FindCrossSectionMax(table);
	}
        if(nullptr == v) { proc->SetCrossSectionType(fEmIncreasing); }
      }
      proc->SetEnergyOfCrossSectionMax(v);
    }
  }
  // protection against double printout
  if(isLocked) { return; }

  // explicitly defined printout by particle name
  if(1 < verb || (0 < verb && (num == "gamma" || num == "e-" || 
			       num == "e+"    || num == "mu+" || 
			       num == "mu-"   || num == "proton"|| 
			       num == "pi+"   || num == "pi-" || 
			       num == "kaon+" || num == "kaon-" || 
			       num == "alpha" || num == "anti_proton" || 
			       num == "GenericIon" ||
			       num == "alpha+" || num == "helium" ||
			       num == "hydrogen"))) { 
    proc->StreamInfo(G4cout, *part);
  }

  if(1 < verb) {
    G4cout << "### G4VEmProcess::BuildPhysicsTable() done for "
           << proc->GetProcessName() << " and particle " << num
           << " baseMat=" << baseMat << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4EmTableUtil::BuildLambdaTable(G4VEmProcess* proc,
                                     const G4ParticleDefinition* part,
                                     G4EmModelManager* modelManager,
			             G4LossTableBuilder* bld,
                                     G4PhysicsTable* theLambdaTable,
                                     G4PhysicsTable* theLambdaTablePrim,
                                     const G4double minKinEnergy,
                                     const G4double minKinEnergyPrim,
                                     const G4double maxKinEnergy,
                                     const G4double scale,
				     const G4int verboseLevel,
                                     const G4bool startFromNull,
                                     const G4bool splineFlag)
{
  if(1 < verboseLevel) {
    G4cout << "G4EmProcess::BuildLambdaTable() for process "
           << proc->GetProcessName() << " and particle "
           << part->GetParticleName() << G4endl;
  }

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t numOfCouples = theCoupleTable->GetTableSize();

  G4PhysicsLogVector* aVector = nullptr;
  G4PhysicsLogVector* aVectorPrim = nullptr;
  G4PhysicsLogVector* bVectorPrim = nullptr;

  G4double emax1 = std::min(maxKinEnergy, minKinEnergyPrim);
    
  for(std::size_t i=0; i<numOfCouples; ++i) {

    if (bld->GetFlag(i)) {
      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = 
        theCoupleTable->GetMaterialCutsCouple((G4int)i);

      // build main table
      if(nullptr != theLambdaTable) {
        delete (*theLambdaTable)[i];

        // if start from zero then change the scale
        G4double emin = minKinEnergy;
        G4bool startNull = false;
        if(startFromNull) {
          G4double e = proc->MinPrimaryEnergy(part, couple->GetMaterial());
          if(e >= emin) {
            emin = e;
            startNull = true;
          }
        }
        G4double emax = emax1;
        if(emax <= emin) { emax = 2*emin; }
        G4int bin = G4lrint(scale*G4Log(emax/emin));
        bin = std::max(bin, 5);
        aVector = new G4PhysicsLogVector(emin, emax, bin, splineFlag);
        modelManager->FillLambdaVector(aVector, couple, startNull);
        if(splineFlag) { aVector->FillSecondDerivatives(); }
        G4PhysicsTableHelper::SetPhysicsVector(theLambdaTable, i, aVector);
      }
      // build high energy table
      if(nullptr != theLambdaTablePrim) {
        delete (*theLambdaTablePrim)[i];

        // start not from zero and always use spline
        if(nullptr == bVectorPrim) {
          G4int bin = G4lrint(scale*G4Log(maxKinEnergy/minKinEnergyPrim));
          bin = std::max(bin, 5);
          aVectorPrim = 
            new G4PhysicsLogVector(minKinEnergyPrim, maxKinEnergy, bin, true);
          bVectorPrim = aVectorPrim;
        } else {
          aVectorPrim = new G4PhysicsLogVector(*bVectorPrim);
        }
        modelManager->FillLambdaVector(aVectorPrim, couple, false, 
                                       fIsCrossSectionPrim);
        aVectorPrim->FillSecondDerivatives();
        G4PhysicsTableHelper::SetPhysicsVector(theLambdaTablePrim, i, 
                                               aVectorPrim);
      }
    }
  }

  if(1 < verboseLevel) {
    G4cout << "Lambda table is built for " << part->GetParticleName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void  G4EmTableUtil::BuildLambdaTable(G4VEnergyLossProcess* proc,
                                     const G4ParticleDefinition* part,
                                     G4EmModelManager* modelManager,
			             G4LossTableBuilder* bld,
                                     G4PhysicsTable* theLambdaTable,
                                     const G4DataVector* theCuts,
                                     const G4double minKinEnergy,
                                     const G4double maxKinEnergy,
                                     const G4double scale,
				     const G4int verboseLevel,
                                     const G4bool splineFlag)
{
  if(1 < verboseLevel) {
    G4cout << "G4EnergyLossProcess::BuildLambdaTable() for process "
           << proc->GetProcessName() << " and particle "
           << part->GetParticleName() << G4endl;
  }

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t numOfCouples = theCoupleTable->GetTableSize();

  G4PhysicsLogVector* aVector = nullptr;
  for(std::size_t i=0; i<numOfCouples; ++i) {
    if (bld->GetFlag(i)) {
      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = 
        theCoupleTable->GetMaterialCutsCouple((G4int)i);

      delete (*theLambdaTable)[i];
      G4bool startNull = true;
      G4double emin = 
        proc->MinPrimaryEnergy(part, couple->GetMaterial(), (*theCuts)[i]);
      if(minKinEnergy > emin) { 
        emin = minKinEnergy; 
        startNull = false;
      }

      G4double emax = maxKinEnergy;
      if(emax <= emin) { emax = 2*emin; }
      G4int bin = G4lrint(scale*G4Log(emax/emin));
      bin = std::max(bin, 5);
      aVector = new G4PhysicsLogVector(emin, emax, bin, splineFlag);
      modelManager->FillLambdaVector(aVector, couple, startNull, fRestricted);
      if(splineFlag) { aVector->FillSecondDerivatives(); }
      G4PhysicsTableHelper::SetPhysicsVector(theLambdaTable, i, aVector);
    }
  }

  if(1 < verboseLevel) {
    G4cout << "Lambda table is built for " << part->GetParticleName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const G4ParticleDefinition*
G4EmTableUtil::CheckIon(G4VEnergyLossProcess* proc,
                        const G4ParticleDefinition* part,
                        const G4ParticleDefinition* partLocal,
                        const G4int verb, G4bool& isIon)
{
  if(1 < verb) {
    G4cout << "G4VEnergyLossProcess::PreparePhysicsTable for "
           << proc->GetProcessName() << " for " << part->GetParticleName() 
           << G4endl;
  }
  const G4ParticleDefinition* particle = partLocal;

  // Are particle defined?
  if(nullptr == particle) { particle = part; }
  if(part->GetParticleType() == "nucleus") {
    G4String pname = part->GetParticleName();
    if(pname != "deuteron" && pname != "triton" &&
       pname != "alpha+"   && pname != "alpha") {

      const G4ParticleDefinition* theGIon = G4GenericIon::GenericIon();
      isIon = true; 

      if(particle != theGIon) {
        G4ProcessManager* pm = theGIon->GetProcessManager();
        G4ProcessVector* v = pm->GetAlongStepProcessVector();
        G4int n = (G4int)v->size();
        for(G4int j=0; j<n; ++j) {
          if((*v)[j] == proc) {
            particle = theGIon;
            break;
          } 
        }
      }
    }
  }
  return particle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4EmTableUtil::UpdateModels(G4VEnergyLossProcess* proc,
			         G4EmModelManager* modelManager,
                                 const G4double maxKinEnergy,
                                 const G4int nModels,
                                 G4int& secID, G4int& biasID,
                                 G4int& mainSec, const G4bool baseMat,
                                 const G4bool isMaster, const G4bool useAGen)
{
  // defined ID of secondary particles
  G4int stype = proc->GetProcessSubType();
  if(stype == fBremsstrahlung) {
    secID = _Bremsstrahlung;
    biasID = _SplitBremsstrahlung;
  } else if(stype == fPairProdByCharged) {
    secID = _PairProduction;
    mainSec = 2;
  }

  // initialisation of models
  for(G4int i=0; i<nModels; ++i) {
    G4VEmModel* mod = modelManager->GetModel(i);
    mod->SetMasterThread(isMaster);
    mod->SetAngularGeneratorFlag(useAGen);
    if(mod->HighEnergyLimit() > maxKinEnergy) {
      mod->SetHighEnergyLimit(maxKinEnergy);
    }
    mod->SetUseBaseMaterials(baseMat);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void 
G4EmTableUtil::BuildLocalElossProcess(G4VEnergyLossProcess* proc,
				      const G4VEnergyLossProcess* masterProc,
				      const G4ParticleDefinition* part,
				      const G4int nModels)
{
  // copy table pointers from master thread
  proc->SetDEDXTable(masterProc->DEDXTable(),fRestricted);
  proc->SetDEDXTable(masterProc->DEDXunRestrictedTable(),fTotal);
  proc->SetDEDXTable(masterProc->IonisationTable(),fIsIonisation);
  proc->SetRangeTableForLoss(masterProc->RangeTableForLoss());
  proc->SetCSDARangeTable(masterProc->CSDARangeTable());
  proc->SetInverseRangeTable(masterProc->InverseRangeTable());
  proc->SetLambdaTable(masterProc->LambdaTable());
  proc->SetCrossSectionType(masterProc->CrossSectionType());
  proc->SetEnergyOfCrossSectionMax(masterProc->EnergyOfCrossSectionMax());
  proc->SetTwoPeaksXS(masterProc->TwoPeaksXS());
  proc->SetIonisation(masterProc->IsIonisationProcess());
  G4bool baseMat = masterProc->UseBaseMaterial();

  // local initialisation of models
  G4bool printing = true;
  for(G4int i=0; i<nModels; ++i) {
    G4VEmModel* mod = proc->GetModelByIndex(i, printing);
    G4VEmModel* mod0= masterProc->GetModelByIndex(i, printing);
    mod->SetUseBaseMaterials(baseMat);
    mod->InitialiseLocal(part, mod0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4EmTableUtil::BuildDEDXTable(G4VEnergyLossProcess* proc,
				   const G4ParticleDefinition* part,
				   G4EmModelManager* modelManager,
				   G4LossTableBuilder* bld,
				   G4PhysicsTable* table,
				   const G4double emin,
				   const G4double emax,
				   const G4int nbins,
				   const G4int verbose,
				   const G4EmTableType tType,
				   const G4bool spline)
{
  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t numOfCouples = theCoupleTable->GetTableSize();

  if(1 < verbose) {
    G4cout << numOfCouples << " couples" << " minKinEnergy(MeV)= " << emin
           << " maxKinEnergy(MeV)= " << emax << " nbins= " << nbins << G4endl;
  }
  G4PhysicsLogVector* aVector = nullptr;
  G4PhysicsLogVector* bVector = nullptr;

  for(std::size_t i=0; i<numOfCouples; ++i) {

    if(1 < verbose) {
      G4cout << "G4VEnergyLossProcess::BuildDEDXVector idx= " << i 
             << "  flagTable=" << table->GetFlag(i) 
             << " flagBuilder=" << bld->GetFlag(i) << G4endl;
    }
    if(bld->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = 
        theCoupleTable->GetMaterialCutsCouple((G4int)i);
      delete (*table)[i];
      if(nullptr != bVector) {
        aVector = new G4PhysicsLogVector(*bVector);
      } else {
        bVector = new G4PhysicsLogVector(emin, emax, nbins, spline);
        aVector = bVector;
      }

      modelManager->FillDEDXVector(aVector, couple, tType);
      if(spline) { aVector->FillSecondDerivatives(); }

      // Insert vector for this material into the table
      G4PhysicsTableHelper::SetPhysicsVector(table, i, aVector);
    }
  }

  if(1 < verbose) {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTable(): table is built for "
           << part->GetParticleName()
           << " and process " << proc->GetProcessName()
           << G4endl;
    if(2 < verbose) G4cout << (*table) << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4EmTableUtil::PrepareMscProcess(G4VMultipleScattering* proc,
				      const G4ParticleDefinition& part,
				      G4EmModelManager* modelManager,
				      G4MscStepLimitType& stepLimit,
                                      G4double& facrange,
				      G4bool& latDisplacement, G4bool& master,
				      G4bool& isIon, G4bool& baseMat)
{
  auto param = G4EmParameters::Instance();
  G4int verb = (master) ? param->Verbose() : param->WorkerVerbose(); 
  proc->SetVerboseLevel(verb);

  if(part.GetPDGMass() > CLHEP::GeV ||
     part.GetParticleName() == "GenericIon") { isIon = true; }

  if(1 < verb) {
    G4cout << "### G4VMultipleScattering::PrepearPhysicsTable() for "
           << proc->GetProcessName()
           << " and particle " << part.GetParticleName()
           << " isIon: " << isIon << " isMaster: " << master
	   << G4endl;
  }

  // initialise process
  proc->InitialiseProcess(&part);

  // heavy particles 
  if(part.GetPDGMass() > CLHEP::MeV) {
    stepLimit = param->MscMuHadStepLimitType(); 
    facrange = param->MscMuHadRangeFactor(); 
    latDisplacement = param->MuHadLateralDisplacement();
  } else {
    stepLimit = param->MscStepLimitType(); 
    facrange = param->MscRangeFactor(); 
    latDisplacement = param->LateralDisplacement();
  }

  // initialisation of models
  auto numberOfModels = modelManager->NumberOfModels();
  for(G4int i=0; i<numberOfModels; ++i) {
    G4VMscModel* msc = proc->GetModelByIndex(i);
    msc->SetIonisation(nullptr, &part);
    msc->SetMasterThread(master);
    msc->SetPolarAngleLimit(param->MscThetaLimit());
    G4double emax = std::min(msc->HighEnergyLimit(),param->MaxKinEnergy());
    msc->SetHighEnergyLimit(emax);
    msc->SetUseBaseMaterials(baseMat);
  }
  modelManager->Initialise(&part, nullptr, verb);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4EmTableUtil::BuildMscProcess(G4VMultipleScattering* proc,
                                    const G4VMultipleScattering* masterProc,
		                    const G4ParticleDefinition& part,
		                    const G4ParticleDefinition* firstPart,
				    G4int nModels, G4bool master)
{
  auto param = G4EmParameters::Instance();
  G4int verb = param->Verbose(); 

  if(!master && firstPart == &part) {
    // initialisation of models
    G4bool baseMat = masterProc->UseBaseMaterial();
    for(G4int i=0; i<nModels; ++i) {
      G4VMscModel* msc = proc->GetModelByIndex(i);
      G4VMscModel* msc0 = masterProc->GetModelByIndex(i);
      msc->SetUseBaseMaterials(baseMat);
      msc->SetCrossSectionTable(msc0->GetCrossSectionTable(), false);
      msc->InitialiseLocal(&part, msc0);
    }
  }
  if(!param->IsPrintLocked()) {
    const G4String& num = part.GetParticleName();

    // explicitly defined printout by particle name
    if(1 < verb || (0 < verb && (num == "e-" || 
		 		 num == "e+"    || num == "mu+" || 
				 num == "mu-"   || num == "proton"|| 
				 num == "pi+"   || num == "pi-" || 
				 num == "kaon+" || num == "kaon-" || 
				 num == "alpha" || num == "anti_proton" || 
				 num == "GenericIon" || num == "alpha+" || 
				 num == "alpha" ))) { 
      proc->StreamInfo(G4cout, part);
    }
  }
  if(1 < verb) {
    G4cout << "### G4VMultipleScattering::BuildPhysicsTable() done for "
	   << proc->GetProcessName()
	   << " and particle " << part.GetParticleName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmTableUtil::StoreMscTable(G4VMultipleScattering* proc,
                                    const G4ParticleDefinition* part,
                                    const G4String& dir,
                                    const G4int nModels, const G4int verb,
		                    const G4bool ascii)
{
  G4bool ok = true;
  for(G4int i=0; i<nModels; ++i) {
    G4VMscModel* msc = proc->GetModelByIndex(i);
    G4PhysicsTable* table = msc->GetCrossSectionTable();
    if (nullptr != table) {
      G4String ss = G4UIcommand::ConvertToString(i);
      G4String name = 
        proc->GetPhysicsTableFileName(part, dir, "LambdaMod"+ss, ascii);
      G4bool yes = table->StorePhysicsTable(name,ascii);

      if ( yes ) {
        if ( verb > 0 ) {
          G4cout << "Physics table are stored for " 
                 << part->GetParticleName()
                 << " and process " << proc->GetProcessName()
                 << " with a name <" << name << "> " << G4endl;
        }
      } else {
        G4cout << "Fail to store Physics Table for " 
               << part->GetParticleName()
               << " and process " << proc->GetProcessName()
               << " in the directory <" << dir
               << "> " << G4endl;
	ok = false;
      }
    }
  }
  return ok;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmTableUtil::StoreTable(G4VProcess* ptr,
                                 const G4ParticleDefinition* part, 
                                 G4PhysicsTable* aTable, 
                                 const G4String& dir,
                                 const G4String& tname,
                                 const G4int verb, const G4bool ascii)
{
  G4bool res = true;
  if (nullptr != aTable) {
    const G4String& name = 
      ptr->GetPhysicsTableFileName(part, dir, tname, ascii);
    if ( aTable->StorePhysicsTable(name, ascii) ) {
      if (1 < verb) G4cout << "Stored: " << name << G4endl;
    } else {
      res = false;
      G4cout << "Fail to store: " << name << G4endl;
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmTableUtil::RetrieveTable(G4VProcess* ptr,
                                    const G4ParticleDefinition* part, 
                                    G4PhysicsTable* aTable, 
                                    const G4String& dir, const G4String& tname,
                                    const G4int verb, const G4bool ascii,
                                    const G4bool spline)
{
  G4bool res = true;
  if (nullptr == aTable) { return res; }
  G4cout << tname << " table for " << part->GetParticleName() 
	 << " will be retrieved " << G4endl;
  const G4String& name = 
    ptr->GetPhysicsTableFileName(part, dir, tname, ascii);
  if(G4PhysicsTableHelper::RetrievePhysicsTable(aTable, name, ascii, spline)) {
    if(spline) {
      for(auto & v : *aTable) {
	if(nullptr != v) { v->FillSecondDerivatives(); }
      }
    }
    if (0 < verb) {
      G4cout << tname << " table for " << part->GetParticleName() 
	     << " is Retrieved from <" << name << ">"
	     << G4endl;
    }
  } else {
    res = false;
    G4cout << "Fail to retrieve: " << tname << " from " << name << " for " 
	   << part->GetParticleName() << G4endl;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

