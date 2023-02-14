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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4GammaGeneralProcess
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 19.07.2018
//
// Modifications:
//
// Class Description:
//

// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4GammaGeneralProcess.hh"
#include "G4VDiscreteProcess.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4ProductionCutsTable.hh"
#include "G4LossTableBuilder.hh"
#include "G4HadronicProcess.hh"
#include "G4LossTableManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4DataVector.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4VParticleChange.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4EmParameters.hh"
#include "G4EmProcessSubType.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4GammaConversionToMuons.hh"
#include "G4Gamma.hh"

#include "G4Log.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmDataHandler* G4GammaGeneralProcess::theHandler = nullptr;
G4bool G4GammaGeneralProcess::theT[nTables] =
  {true,false,true,true,true,false,true,true,true,
   true,true,true,true,true,true};
G4String G4GammaGeneralProcess::nameT[nTables] =
  {"0","1","2","3","4","5","6","7","8",
   "9","10","11","12","13","14"};

G4GammaGeneralProcess::G4GammaGeneralProcess(const G4String& pname):
  G4VEmProcess(pname, fElectromagnetic),
  minPEEnergy(150*CLHEP::keV),
  minEEEnergy(2*CLHEP::electron_mass_c2),
  minMMEnergy(100*CLHEP::MeV)
{
  SetVerboseLevel(1);
  SetParticle(G4Gamma::Gamma());
  SetProcessSubType(fGammaGeneralProcess);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4GammaGeneralProcess::~G4GammaGeneralProcess()
{
  if(isTheMaster) {
    delete theHandler;
    theHandler = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4GammaGeneralProcess::IsApplicable(const G4ParticleDefinition&)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4GammaGeneralProcess::AddEmProcess(G4VEmProcess* ptr)
{
  if(nullptr == ptr) { return; }
  G4int stype = ptr->GetProcessSubType();
  if(stype == fRayleigh)                 { theRayleigh = ptr; }
  else if(stype == fPhotoElectricEffect) { thePhotoElectric = ptr; }
  else if(stype == fComptonScattering)   { theCompton = ptr; }
  else if(stype == fGammaConversion)     { theConversionEE = ptr; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4GammaGeneralProcess::AddMMProcess(G4GammaConversionToMuons* ptr)
{
  theConversionMM = ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4GammaGeneralProcess::AddHadProcess(G4HadronicProcess* ptr)
{
  theGammaNuclear = ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4GammaGeneralProcess::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  SetParticle(&part);
  preStepLambda = 0.0;
  idxEnergy = 0;
  currentCouple = nullptr;

  G4EmParameters* param = G4EmParameters::Instance();
  G4LossTableManager* man = G4LossTableManager::Instance();

  isTheMaster = man->IsMaster();
  if(isTheMaster) { SetVerboseLevel(param->Verbose()); }
  else { SetVerboseLevel(param->WorkerVerbose()); }

  G4LossTableBuilder* bld = man->GetTableBuilder();
  baseMat = bld->GetBaseMaterialFlag();

  if(1 < verboseLevel) {
    G4cout << "G4GammaGeneralProcess::PreparePhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
	   << " isMaster: " << isTheMaster << G4endl;
  }

  // 3 sub-processes must be always defined
  if(thePhotoElectric == nullptr || theCompton == nullptr ||
     theConversionEE == nullptr) {
    G4ExceptionDescription ed;
    ed << "### G4GeneralGammaProcess is initialized incorrectly"
       << "\n Photoelectric: " << thePhotoElectric
       << "\n Compton: " << theCompton
       << "\n Conversion: " << theConversionEE;
    G4Exception("G4GeneralGammaProcess","em0004",
		FatalException, ed,"");
  }

  if(thePhotoElectric) { thePhotoElectric->PreparePhysicsTable(part); }
  if(theCompton)       { theCompton->PreparePhysicsTable(part); }
  if(theConversionEE)  { theConversionEE->PreparePhysicsTable(part); }
  if(theRayleigh)      { theRayleigh->PreparePhysicsTable(part); }
  if(theGammaNuclear)  { theGammaNuclear->PreparePhysicsTable(part); }
  if(theConversionMM)  { theConversionMM->PreparePhysicsTable(part); }

  InitialiseProcess(&part);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4GammaGeneralProcess::InitialiseProcess(const G4ParticleDefinition*)
{
  if(isTheMaster) {

    G4EmParameters* param = G4EmParameters::Instance();
    G4LossTableManager* man = G4LossTableManager::Instance();

    // tables are created and its size is defined only once
    if(nullptr == theHandler) {
      theHandler = new G4EmDataHandler(nTables);
      if(theRayleigh) { theT[1] = true; }

      theHandler->SetMasterProcess(thePhotoElectric);
      theHandler->SetMasterProcess(theCompton);
      theHandler->SetMasterProcess(theConversionEE);
      theHandler->SetMasterProcess(theRayleigh);
    }
    auto bld = man->GetTableBuilder();

    const G4ProductionCutsTable* theCoupleTable=
      G4ProductionCutsTable::GetProductionCutsTable();
    std::size_t numOfCouples = theCoupleTable->GetTableSize();

    G4double mine = param->MinKinEnergy();
    G4double maxe = param->MaxKinEnergy();
    G4int nd = param->NumberOfBinsPerDecade();
    std::size_t nbin1 = std::max(5, nd*G4lrint(std::log10(minPEEnergy/mine)));
    std::size_t nbin2 = std::max(5, nd*G4lrint(std::log10(maxe/minMMEnergy)));

    G4PhysicsVector* vec = nullptr;
    G4PhysicsLogVector aVector(mine,minPEEnergy,nbin1,true);
    G4PhysicsLogVector bVector(minPEEnergy,minEEEnergy,nLowE,false);
    G4PhysicsLogVector cVector(minEEEnergy,minMMEnergy,nHighE,false);
    G4PhysicsLogVector dVector(minMMEnergy,maxe,nbin2,true);

    for(std::size_t i=0; i<nTables; ++i) {
      if(!theT[i]) { continue; }
      //G4cout << "## PreparePhysTable " << i << "." << G4endl;
      G4PhysicsTable* table = theHandler->MakeTable(i);
      //G4cout << "   make table " << table << G4endl;
      for(std::size_t j=0; j<numOfCouples; ++j) {
	vec = (*table)[j];
	if (bld->GetFlag(j) && nullptr == vec) {
	  //G4cout <<"InitialiseProcess iTable="<<i<<" jCouple="<< j <<" make new vector"<< G4endl;
	  if(i<=1) {
	    vec = new G4PhysicsVector(aVector);
	  } else if(i<=5) {
	    vec = new G4PhysicsVector(bVector);
	  } else if(i<=9) {
	    vec = new G4PhysicsVector(cVector);
	  } else {
	    vec = new G4PhysicsVector(dVector);
	  }
	  G4PhysicsTableHelper::SetPhysicsVector(table, j, vec);
	}
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4GammaGeneralProcess::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  if(1 < verboseLevel) {
    G4cout << "### G4VEmProcess::BuildPhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << G4endl;
  }
  if(!isTheMaster) {
    thePhotoElectric->SetEmMasterProcess(theHandler->GetMasterProcess(0));
    baseMat = theHandler->GetMasterProcess(0)->UseBaseMaterial();
  }
  thePhotoElectric->BuildPhysicsTable(part);

  if(!isTheMaster) {
    theCompton->SetEmMasterProcess(theHandler->GetMasterProcess(1));
  }
  theCompton->BuildPhysicsTable(part);

  if(!isTheMaster) {
    theConversionEE->SetEmMasterProcess(theHandler->GetMasterProcess(2));
  }
  theConversionEE->BuildPhysicsTable(part);

  if(theRayleigh != nullptr) {
    if(!isTheMaster) {
      theRayleigh->SetEmMasterProcess(theHandler->GetMasterProcess(3));
    }
    theRayleigh->BuildPhysicsTable(part);
  }
  if(theGammaNuclear != nullptr)  { theGammaNuclear->BuildPhysicsTable(part); }
  if(theConversionMM != nullptr)  { theConversionMM->BuildPhysicsTable(part); }

  if(isTheMaster) {
    const G4ProductionCutsTable* theCoupleTable=
      G4ProductionCutsTable::GetProductionCutsTable();
    G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

    G4LossTableBuilder* bld = G4LossTableManager::Instance()->GetTableBuilder();
    const std::vector<G4PhysicsTable*>& tables = theHandler->GetTables();

    G4CrossSectionDataStore* gn = (nullptr != theGammaNuclear)
      ? theGammaNuclear->GetCrossSectionDataStore() : nullptr;
    G4DynamicParticle* dynParticle =
      new G4DynamicParticle(G4Gamma::Gamma(),G4ThreeVector(1,0,0),1.0);

    G4double sigComp(0.), sigPE(0.), sigConv(0.), sigR(0.),
      sigN(0.), sigM(0.), val(0.);

    for(G4int i=0; i<numOfCouples; ++i) {

      if (bld->GetFlag(i)) {
        G4int idx = (!baseMat) ? i : DensityIndex(i);
	const G4MaterialCutsCouple* couple =
	  theCoupleTable->GetMaterialCutsCouple(i);
	const G4Material* material = couple->GetMaterial();

	// energy interval 0
        std::size_t nn = (*(tables[0]))[idx]->GetVectorLength();
	if(1 < verboseLevel) {
          G4cout << "======= Zone 0 ======= N= " << nn
		 << " for " << material->GetName() << G4endl;
	}
        for(std::size_t j=0; j<nn; ++j) {
          G4double e = (*(tables[0]))[idx]->Energy(j);
          G4double loge = G4Log(e);
          sigComp = theCompton->GetLambda(e, couple, loge);
          sigR = (nullptr != theRayleigh) ?
            theRayleigh->GetLambda(e, couple, loge) : 0.0;
          G4double sum = sigComp + sigR;
	  if(1 < verboseLevel) {
	    G4cout << j << ". E= " << e << " xs= " << sum
		   << " compt= " << sigComp << " Rayl= " << sigR << G4endl;
	  }
          (*(tables[0]))[idx]->PutValue(j, sum);
	  if(theT[1]) {
            val = sigR/sum;
	    (*(tables[1]))[idx]->PutValue(j, val);
	  }
	}

	// energy interval 1
        nn = (*(tables[2]))[idx]->GetVectorLength();
	if(1 < verboseLevel) {
	  G4cout << "======= Zone 1 ======= N= " << nn << G4endl;
	}
        for(std::size_t j=0; j<nn; ++j) {
          G4double e = (*(tables[2]))[idx]->Energy(j);
          G4double loge = G4Log(e);
          sigComp = theCompton->GetLambda(e, couple, loge);
          sigR = (nullptr != theRayleigh) ?
            theRayleigh->GetLambda(e, couple, loge) : 0.0;
          sigPE = thePhotoElectric->GetLambda(e, couple, loge);
          G4double sum = sigComp + sigR + sigPE;
	  if(1 < verboseLevel) {
	    G4cout << j << ". E= " << e << " xs= " << sum
		   << " compt= " << sigComp << " conv= " << sigConv
		   << " PE= " << sigPE << " Rayl= " << sigR
		   << " GN= " << sigN << G4endl;
	  }
          (*(tables[2]))[idx]->PutValue(j, sum);

          val = sigPE/sum;
	  (*(tables[3]))[idx]->PutValue(j, val);

	  val = (sigR > 0.0) ? (sigComp + sigPE)/sum : 1.0;
	  (*(tables[4]))[idx]->PutValue(j, val);
	}

	// energy interval 2
        nn = (*(tables[6]))[idx]->GetVectorLength();
	if(1 < verboseLevel) {
	  G4cout << "======= Zone 2 ======= N= " << nn << G4endl;
	}
        for(std::size_t j=0; j<nn; ++j) {
          G4double e = (*(tables[6]))[idx]->Energy(j);
          G4double loge = G4Log(e);
          sigComp = theCompton->GetLambda(e, couple, loge);
          sigConv = theConversionEE->GetLambda(e, couple, loge);
          sigPE = thePhotoElectric->GetLambda(e, couple, loge);
          sigN = 0.0;
          if(nullptr != gn) {
	    dynParticle->SetKineticEnergy(e);
	    sigN = gn->ComputeCrossSection(dynParticle, material);
	  }
          G4double sum = sigComp + sigConv + sigPE + sigN;
	  if(1 < verboseLevel) {
	    G4cout << j << ". E= " << e << " xs= " << sum
		   << " compt= " << sigComp << " conv= " << sigConv
		   << " PE= " << sigPE
		   << " GN= " << sigN << G4endl;
	  }
          (*(tables[6]))[idx]->PutValue(j, sum);

          val = sigConv/sum;
	  (*(tables[7]))[idx]->PutValue(j, val);

          val = (sigConv + sigComp)/sum;
	  (*(tables[8]))[idx]->PutValue(j, val);

	  val = (sigN > 0.0) ? (sigConv + sigComp + sigPE)/sum : 1.0;
	  (*(tables[9]))[idx]->PutValue(j, val);
	}

	// energy interval 3
        nn = (*(tables[10]))[idx]->GetVectorLength();
	if(1 < verboseLevel) {
	  G4cout << "======= Zone 3 ======= N= " << nn
		 << " for " << material->GetName() << G4endl;
	}
        for(std::size_t j=0; j<nn; ++j) {
          G4double e = (*(tables[10]))[idx]->Energy(j);
          G4double loge = G4Log(e);
          sigComp = theCompton->GetLambda(e, couple, loge);
          sigConv = theConversionEE->GetLambda(e, couple, loge);
          sigPE = thePhotoElectric->GetLambda(e, couple, loge);
          sigN = 0.0;
          if(nullptr != gn) {
	    dynParticle->SetKineticEnergy(e);
	    sigN = gn->ComputeCrossSection(dynParticle, material);
	  }
          sigM = 0.0;
          if(nullptr != theConversionMM) {
	    val = theConversionMM->ComputeMeanFreePath(e, material);
	    sigM = (val < DBL_MAX) ? 1./val : 0.0;
	  }
          G4double sum = sigComp + sigConv + sigPE + sigN + sigM;
	  if(1 < verboseLevel) {
	    G4cout << j << ". E= " << e << " xs= " << sum
		   << " compt= " << sigComp << " conv= " << sigConv
		   << " PE= " << sigPE
		   << " GN= " << sigN << G4endl;
	  }
          (*(tables[10]))[idx]->PutValue(j, sum);

          val = (sigComp + sigPE + sigN + sigM)/sum;
	  (*(tables[11]))[idx]->PutValue(j, val);

          val = (sigPE + sigN + sigM)/sum;
	  (*(tables[12]))[idx]->PutValue(j, val);

	  val = (sigN + sigM)/sum;
	  (*(tables[13]))[idx]->PutValue(j, val);

	  val = sigN/sum;
	  (*(tables[14]))[idx]->PutValue(j, val);
	}
	for(std::size_t k=0; k<nTables; ++k) {
	  if(theT[k] && (k <= 1 || k >= 10)) {
	    //G4cout <<"BuildPhysTable spline iTable="<<k<<" jCouple="<< idx << G4endl;
	    (*(tables[k]))[idx]->FillSecondDerivatives();
	  }
	}
      }
    }
    delete dynParticle;
  }

  if(1 < verboseLevel) {
    G4cout << "### G4VEmProcess::BuildPhysicsTable() done for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4GammaGeneralProcess::StartTracking(G4Track*)
{
  theNumberOfInteractionLengthLeft = -1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4GammaGeneralProcess::PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition)
{
  *condition = NotForced;
  G4double x = DBL_MAX;

  G4double energy = track.GetKineticEnergy();
  const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();

  // compute mean free path
  G4bool recompute = false;
  if(couple != currentCouple) {
    currentCouple = couple;
    basedCoupleIndex = currentCoupleIndex = couple->GetIndex();
    currentMaterial = couple->GetMaterial();
    factor = 1.0;
    if(baseMat) {
      basedCoupleIndex = DensityIndex((G4int)currentCoupleIndex);
      factor = DensityFactor((G4int)currentCoupleIndex);
    }
    recompute = true;
  }
  if(energy != preStepKinEnergy) {
    preStepKinEnergy = energy;
    preStepLogE = track.GetDynamicParticle()->GetLogKineticEnergy();
    recompute = true;
  }
  if(recompute) {
    preStepLambda = TotalCrossSectionPerVolume();

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
  /*
  G4cout << "PostStepGetPhysicalInteractionLength: e= " << energy
	 << " idxe= " << idxEnergy << "  xs= " << preStepLambda
	 << " x= " << x << G4endl;
  */
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4GammaGeneralProcess::TotalCrossSectionPerVolume()
{
  G4double cross = 0.0;
  /*
  G4cout << "#Total: " << preStepKinEnergy << " " << minPEEnergy << " "
         << minEEEnergy << " " << minMMEnergy<< G4endl;
  G4cout << " idxE= " << idxEnergy
	 << " idxC= " << currentCoupleIndex << G4endl;
  */
  if(preStepKinEnergy < minPEEnergy) {
    cross = ComputeGeneralLambda(0, 0);
    //G4cout << "XS1: " << cross << G4endl;
    peLambda = thePhotoElectric->GetLambda(preStepKinEnergy, currentCouple, preStepLogE);
    cross += peLambda;
    //G4cout << "XS2: " << peLambda << G4endl;

  } else if(preStepKinEnergy < minEEEnergy) {
    cross = ComputeGeneralLambda(1, 2);
    //G4cout << "XS3: " << cross << G4endl;

  } else if(preStepKinEnergy < minMMEnergy) {
    cross = ComputeGeneralLambda(2, 6);
    //G4cout << "XS4: " << cross << G4endl;

  } else {
    cross = ComputeGeneralLambda(3, 10);
    //G4cout << "XS5: " << cross << G4endl;
  }
  /*  
  G4cout << "xs= " << cross << " idxE= " << idxEnergy
	 << " idxC= " << currentCoupleIndex
	 << " E= " << preStepKinEnergy << G4endl;
  */
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4GammaGeneralProcess::PostStepDoIt(const G4Track& track,
                                                       const G4Step& step)
{
  // In all cases clear number of interaction lengths
  theNumberOfInteractionLengthLeft = -1.0;
  selectedProc = nullptr;
  G4double q = G4UniformRand();
  /*
  G4cout << "PostStep: preStepLambda= " << preStepLambda
         << " PE= " << peLambda << " q= " << q << " idxE= " << idxEnergy
         << G4endl;
  */
  switch (idxEnergy) {
  case 0:
    if(preStepLambda*q <= peLambda) {
      SelectEmProcess(step, thePhotoElectric);
    } else {
      if(theT[1] && preStepLambda*q < (preStepLambda - peLambda)*GetProbability(1) + peLambda) {
	SelectEmProcess(step, theRayleigh);
      } else {
	SelectEmProcess(step, theCompton);
      }
    }
    break;

  case 1:
    if(q <= GetProbability(3)) {
      SelectEmProcess(step, thePhotoElectric);
    } else if(q <= GetProbability(4)) {
      SelectEmProcess(step, theCompton);
    } else if(theRayleigh) {
      SelectEmProcess(step, theRayleigh);
    } else {
      SelectEmProcess(step, thePhotoElectric);
    }
    break;

  case 2:
    if(q <= GetProbability(7)) {
      SelectEmProcess(step, theConversionEE);
    } else if(q <= GetProbability(8)) {
      SelectEmProcess(step, theCompton);
    } else if(q <= GetProbability(9)) {
      SelectEmProcess(step, thePhotoElectric);
    } else if(theGammaNuclear) {
      SelectHadProcess(track, step, theGammaNuclear);
    } else {
      SelectEmProcess(step, theConversionEE);
    }
    break;

  case 3:
    if(q + GetProbability(11) <= 1.0) {
      SelectEmProcess(step, theConversionEE);
    } else if(q + GetProbability(12) <= 1.0) {
      SelectEmProcess(step, theCompton);
    } else if(q + GetProbability(13) <= 1.0) {
      SelectEmProcess(step, thePhotoElectric);
    } else if(theGammaNuclear && q + GetProbability(14) <= 1.0) {
      SelectHadProcess(track, step, theGammaNuclear);
    } else if(theConversionMM) {
      SelectedProcess(step, theConversionMM);
    } else {
      SelectEmProcess(step, theConversionEE);
    }
    break;
  }
  // sample secondaries
  if(selectedProc != nullptr) {
    return selectedProc->PostStepDoIt(track, step);
  }
  // no interaction - exception case
  fParticleChange.InitializeForPostStep(track);
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4GammaGeneralProcess::SelectHadProcess(const G4Track& track,
            const G4Step& step, G4HadronicProcess* proc)
{
  SelectedProcess(step, proc);
  proc->GetCrossSectionDataStore()->ComputeCrossSection(track.GetDynamicParticle(),
                                                        currentMaterial);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4GammaGeneralProcess::StorePhysicsTable(const G4ParticleDefinition* part,
                                              const G4String& directory,
                                              G4bool ascii)
{
  G4bool yes = true;
  if(!isTheMaster) { return yes; }
  if(!thePhotoElectric->StorePhysicsTable(part, directory, ascii))
    { yes = false; }
  if(!theCompton->StorePhysicsTable(part, directory, ascii))
    { yes = false; }
  if(!theConversionEE->StorePhysicsTable(part, directory, ascii))
    { yes = false; }
  if(theRayleigh != nullptr &&
     !theRayleigh->StorePhysicsTable(part, directory, ascii))
    { yes = false; }

  for(std::size_t i=0; i<nTables; ++i) {
    if(theT[i]) {
      G4String nam = (0==i || 2==i || 6==i || 10==i)
	? "LambdaGeneral" + nameT[i] : "ProbGeneral" + nameT[i];
      G4String fnam =  GetPhysicsTableFileName(part,directory,nam,ascii);
      if(!theHandler->StorePhysicsTable(i, part, fnam, ascii)) { yes = false; }
    }
  }
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool
G4GammaGeneralProcess::RetrievePhysicsTable(const G4ParticleDefinition* part,
                                            const G4String& directory,
                                            G4bool ascii)
{
  if(1 < verboseLevel) {
    G4cout << "G4GammaGeneralProcess::RetrievePhysicsTable() for "
           << part->GetParticleName() << " and process "
           << GetProcessName() << G4endl;
  }
  G4bool yes = true;
  if(!thePhotoElectric->RetrievePhysicsTable(part, directory, ascii))
    { yes = false; }
  if(!theCompton->RetrievePhysicsTable(part, directory, ascii))
    { yes = false; }
  if(!theConversionEE->RetrievePhysicsTable(part, directory, ascii))
    { yes = false; }
  if(theRayleigh != nullptr &&
     !theRayleigh->RetrievePhysicsTable(part, directory, ascii))
    { yes = false; }

  for(std::size_t i=0; i<nTables; ++i) {
    if(theT[i]) {
      G4String nam = (0==i || 2==i || 6==i || 10==i)
	? "LambdaGeneral" + nameT[i] : "ProbGeneral" + nameT[i];
      G4String fnam = GetPhysicsTableFileName(part,directory,nam,ascii);
      G4bool spline = (i <= 1 || i >= 10); 
      if(!theHandler->RetrievePhysicsTable(i, part, fnam, ascii, spline))
	{ yes = false; }
    }
  }
  if(yes) {
  }
  return yes;
}

//....Ooooo0ooooo ........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4GammaGeneralProcess::GetMeanFreePath(const G4Track& track,
                                              G4double,
                                              G4ForceCondition* condition)
{
  *condition = NotForced;
  return MeanFreePath(track);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4GammaGeneralProcess::ProcessDescription(std::ostream& out) const
{
  thePhotoElectric->ProcessDescription(out);
  theCompton->ProcessDescription(out);
  theConversionEE->ProcessDescription(out);
  if(theRayleigh)      { theRayleigh->ProcessDescription(out); }
  if(theGammaNuclear)  { theGammaNuclear->ProcessDescription(out); }
  if(theConversionMM)  { theConversionMM->ProcessDescription(out); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4String& G4GammaGeneralProcess::GetSubProcessName() const
{
  return (selectedProc) ? selectedProc->GetProcessName()
    : G4VProcess::GetProcessName();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4GammaGeneralProcess::GetSubProcessSubType() const
{
  return (selectedProc) ? selectedProc->GetProcessSubType()
    : fGammaGeneralProcess;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmProcess* G4GammaGeneralProcess::GetEmProcess(const G4String& name)
{
  G4VEmProcess* proc = nullptr;
  if(name == thePhotoElectric->GetProcessName()) {
    proc = thePhotoElectric;
  } else if(name == theCompton->GetProcessName()) {
    proc = theCompton;
  } else if(name == theConversionEE->GetProcessName()) {
    proc = theConversionEE;
  } else if(theRayleigh != nullptr && name == theRayleigh->GetProcessName()) {
    proc = theRayleigh;
  }
  return proc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4VProcess* G4GammaGeneralProcess::GetCreatorProcess() const
{
  return selectedProc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
