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

#include "G4Log.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmDataHandler* G4GammaGeneralProcess::theHandler = nullptr;
G4bool G4GammaGeneralProcess::theT[nTables] = 
  {true,false,true,true,false,false,true,true,true,
   false,true,true,true,false,false};
G4String G4GammaGeneralProcess::nameT[nTables] = 
  {"0","1","2","3","4","5","6","7","8",
   "9","10","11","12","13","14"};

G4GammaGeneralProcess::G4GammaGeneralProcess():
  G4VEmProcess("GammaGeneralProc", fElectromagnetic),
  minPEEnergy(50*CLHEP::keV),
  minEEEnergy(2*CLHEP::electron_mass_c2),
  minMMEnergy(100*CLHEP::MeV),
  peLambda(0.0),
  nLowE(100),
  nHighE(100),
  splineFlag(false)
{
  thePhotoElectric = theCompton = theConversionEE = theRayleigh = nullptr;
  theGammaNuclear = nullptr;
  theConversionMM = nullptr;
  selectedProc = nullptr;

  idxEnergy = idx0 = idx1 = idx2 = idx3 = 0;

  SetVerboseLevel(1);
  SetParticle(theGamma);
  SetProcessSubType(fGammaGeneralProcess);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4GammaGeneralProcess::~G4GammaGeneralProcess()
{
  //std::cout << "G4GammaGeneralProcess::~G4GammaGeneralProcess " << isTheMaster
  //	    << "  " << theHandler << G4endl;
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
  if(!ptr) { return; }
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
  if(1 < verboseLevel) {
    G4cout << "G4GammaGeneralProcess::PreparePhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << G4endl;
  }
  SetParticle(&part);
  currentCouple = nullptr;
  currentMaterial = nullptr;
  preStepLambda = 0.0;
  idxEnergy = idx0 = idx1 = idx2 = idx3 = 0;

  isTheMaster = lManager->IsMaster(); 
  if(isTheMaster) { SetVerboseLevel(theParameters->Verbose()); }
  else { SetVerboseLevel(theParameters->WorkerVerbose()); }

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

    // tables are created and its size is defined only once
    if(!theHandler) { 
      theHandler = new G4EmDataHandler(nTables);  
      if(theRayleigh) { theT[1] = theT[4] = true; } 
      if(theGammaNuclear) { theT[4] = theT[5] = theT[9] = theT[13] = true; } 
      if(theConversionMM) { theT[14] = true; } 

      theHandler->SetMasterProcess(thePhotoElectric);
      theHandler->SetMasterProcess(theCompton);
      theHandler->SetMasterProcess(theConversionEE);
      theHandler->SetMasterProcess(theRayleigh);
     
      const G4ProductionCutsTable* theCoupleTable=
	G4ProductionCutsTable::GetProductionCutsTable();
      size_t numOfCouples = theCoupleTable->GetTableSize();

      G4double mine = theParameters->MinKinEnergy();
      G4double maxe = theParameters->MaxKinEnergy();
      G4int nd = theParameters->NumberOfBinsPerDecade();
      size_t nbin1 = 
	std::max(5, nd*G4lrint(std::log10(minPEEnergy/mine)));
      size_t nbin2 = 
	std::max(5, nd*G4lrint(std::log10(maxe/minMMEnergy)));

      G4PhysicsVector* vec = nullptr;
      G4PhysicsLogVector aVector(mine,minPEEnergy,nbin1);
      G4PhysicsLinearVector bVector(minPEEnergy,minEEEnergy,nLowE);
      G4PhysicsLinearVector cVector(minEEEnergy,minMMEnergy,nHighE);
      G4PhysicsLogVector dVector(minMMEnergy,maxe,nbin2);
      if(splineFlag) { 
        aVector.SetSpline(splineFlag);
        bVector.SetSpline(splineFlag);
        cVector.SetSpline(splineFlag);
        dVector.SetSpline(splineFlag);
      }

      for(size_t i=0; i<nTables; ++i) { 
	if(theT[i]) { 
	  G4PhysicsTable* table = theHandler->MakeTable(i);
          for(size_t j=0; j<numOfCouples; ++j) {
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
  if(thePhotoElectric) { 
    if(!isTheMaster) { 
      thePhotoElectric->SetEmMasterProcess(theHandler->GetMasterProcess(0)); 
    }
    thePhotoElectric->BuildPhysicsTable(part); 
  }
  if(theCompton) { 
    if(!isTheMaster) { 
      theCompton->SetEmMasterProcess(theHandler->GetMasterProcess(1)); 
    }
    theCompton->BuildPhysicsTable(part); 
  }
  if(theConversionEE) { 
    if(!isTheMaster) { 
      theConversionEE->SetEmMasterProcess(theHandler->GetMasterProcess(2)); 
    }
    theConversionEE->BuildPhysicsTable(part); 
  }
  if(theRayleigh) { 
    if(!isTheMaster) { 
      theRayleigh->SetEmMasterProcess(theHandler->GetMasterProcess(3)); 
    }
    theRayleigh->BuildPhysicsTable(part); 
  }
  if(theGammaNuclear)  { theGammaNuclear->BuildPhysicsTable(part); }
  if(theConversionMM)  { theConversionMM->BuildPhysicsTable(part); }

  if(isTheMaster) {
    const G4ProductionCutsTable* theCoupleTable=
      G4ProductionCutsTable::GetProductionCutsTable();
    size_t numOfCouples = theCoupleTable->GetTableSize();

    G4LossTableBuilder* bld = lManager->GetTableBuilder();
    const std::vector<G4PhysicsTable*>& tables = theHandler->GetTables();

    G4CrossSectionDataStore* gn = (theGammaNuclear) 
      ? theGammaNuclear->GetCrossSectionDataStore() : nullptr;
    G4DynamicParticle* dynParticle = 
      new G4DynamicParticle(G4Gamma::Gamma(),G4ThreeVector(1,0,0),1.0);

    G4double sigComp(0.), sigPE(0.), sigConv(0.), sigR(0.),
      sigN(0.), sigM(0.), val(0.);

    for(size_t i=0; i<numOfCouples; ++i) {

      if (bld->GetFlag(i)) {

	const G4MaterialCutsCouple* couple = 
	  theCoupleTable->GetMaterialCutsCouple(i);
	const G4Material* material = couple->GetMaterial();

	// energy interval 0
        size_t nn = (*(tables[0]))[i]->GetVectorLength();
	if(1 < verboseLevel) {
	  G4cout << "======= Zone 0 ======= N= " << nn << G4endl; 
	}
        for(size_t j=0; j<nn; ++j) {
          G4double e = (*(tables[0]))[i]->Energy(j);
          sigComp = (theCompton) ? theCompton->GetLambda(e, couple) : 0.0;
          sigR = (theRayleigh) ? theRayleigh->GetLambda(e, couple) : 0.0;
          G4double sum = sigComp + sigR;
	  if(1 < verboseLevel) {
	    G4cout << j << ". E= " << e << " xs= " << sum 
		   << " compt= " << sigComp << " Rayl= " << sigR << G4endl; 
	  }
          (*(tables[0]))[i]->PutValue(j, sum);
	  if(theT[1]) {
            val = (sum > 0.0) ? sigComp/sum : 0.0;
	    (*(tables[1]))[i]->PutValue(j, val);
	  } 
	}

	// energy interval 1
        nn = (*(tables[2]))[i]->GetVectorLength();
	if(1 < verboseLevel) {
	  G4cout << "======= Zone 1 ======= N= " << nn << G4endl; 
	}
        for(size_t j=0; j<nn; ++j) {
          G4double e = (*(tables[2]))[i]->Energy(j);
          sigComp = (theCompton) ? theCompton->GetLambda(e, couple) : 0.0;
          sigR = (theRayleigh) ? theRayleigh->GetLambda(e, couple) : 0.0;
          sigPE = (thePhotoElectric) 
	    ? thePhotoElectric->GetLambda(e, couple) : 0.0;
          sigN = 0.0;
          if(gn) {
	    dynParticle->SetKineticEnergy(e);
	    sigN = gn->ComputeCrossSection(dynParticle, material);
	  }
          G4double sum = sigComp + sigR + sigPE + sigN;
	  if(1 < verboseLevel) {
	    G4cout << j << ". E= " << e << " xs= " << sum 
		   << " compt= " << sigComp << " conv= " << sigConv 
		   << " PE= " << sigPE << " Rayl= " << sigR
		   << " GN= " << sigN << G4endl; 
	  }
          (*(tables[2]))[i]->PutValue(j, sum);
          val = (sum > 0.0) ? sigPE/sum : 0.0;
	  (*(tables[3]))[i]->PutValue(j, val);
	  if(theT[4]) {
            val = (sum > 0.0) ? (sigComp + sigPE)/sum : 0.0;
	    (*(tables[4]))[i]->PutValue(j, val);
	  } 
	  if(theT[5]) {
            val = (sum > 0.0) ? (sigComp + sigPE + sigR)/sum : 0.0;
	    (*(tables[5]))[i]->PutValue(j, val);
	  } 
	}

	// energy interval 2
        nn = (*(tables[6]))[i]->GetVectorLength();
	if(1 < verboseLevel) {
	  G4cout << "======= Zone 2 ======= N= " << nn << G4endl; 
	}
        for(size_t j=0; j<nn; ++j) {
          G4double e = (*(tables[6]))[i]->Energy(j);
          sigComp = (theCompton) ? theCompton->GetLambda(e, couple) : 0.0;
          sigConv = (theConversionEE) 
	    ? theConversionEE->GetLambda(e, couple) : 0.0;
          sigPE = (thePhotoElectric) 
	    ? thePhotoElectric->GetLambda(e, couple) : 0.0;
          sigN = 0.0;
          if(gn) {
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
          (*(tables[6]))[i]->PutValue(j, sum);
          val = (sum > 0.0) ? sigConv/sum : 0.0;
	  (*(tables[7]))[i]->PutValue(j, val);
          val = (sum > 0.0) ? (sigConv + sigComp)/sum : 0.0;
	  (*(tables[8]))[i]->PutValue(j, val);
	  if(theT[9]) {
            val = (sum > 0.0) ? (sigConv + sigComp + sigPE)/sum : 0.0;
	    (*(tables[9]))[i]->PutValue(j, val);
	  } 
	}

	// energy interval 3
        nn = (*(tables[10]))[i]->GetVectorLength();
	if(1 < verboseLevel) {
	  G4cout << "======= Zone 3 ======= N= " << nn << G4endl; 
	}
        for(size_t j=0; j<nn; ++j) {
          G4double e = (*(tables[10]))[i]->Energy(j);
          sigComp = (theCompton) ? theCompton->GetLambda(e, couple) : 0.0;
          sigConv = (theConversionEE) 
	    ? theConversionEE->GetLambda(e, couple) : 0.0;
          sigPE = (thePhotoElectric) 
	    ? thePhotoElectric->GetLambda(e, couple) : 0.0;
          sigN = 0.0;
          if(gn) {
	    dynParticle->SetKineticEnergy(e);
	    sigN = gn->ComputeCrossSection(dynParticle, material);
	  }
          sigM = 0.0;
          if(theConversionMM) {
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
          (*(tables[10]))[i]->PutValue(j, sum);
          val = (sum > 0.0) ? 1.0 - sigConv/sum : 1.0;
	  (*(tables[11]))[i]->PutValue(j, val);
          val = (sum > 0.0) ? 1.0 - (sigConv + sigComp)/sum : 1.0;
	  (*(tables[12]))[i]->PutValue(j, val);
	  if(theT[13]) {
            val = (sum > 0.0) ? 1.0 - (sigConv + sigComp + sigPE)/sum : 1.0;
	    (*(tables[13]))[i]->PutValue(j, val);
	  } 
	  if(theT[14]) {
            val = (sum > 0.0) 
	      ? 1.0 - (sigConv + sigComp + sigPE + sigN)/sum : 1.0;
	    (*(tables[14]))[i]->PutValue(j, val);
	  } 
	}
	for(size_t k=0; k<nTables; ++k) {
	  if(theT[k] && splineFlag) { 
	    (*(tables[k]))[i]->FillSecondDerivatives(); 
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
  currentMaterial = nullptr;
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
  const G4Material* mat = couple->GetMaterial();

  // compute mean free path
  if(mat != currentMaterial || energy != preStepKinEnergy) {
    currentMaterial = mat;
    preStepKinEnergy = energy;
    preStepLambda = TotalCrossSectionPerVolume(energy, couple);

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

G4double G4GammaGeneralProcess::TotalCrossSectionPerVolume(G4double energy,
                                      const G4MaterialCutsCouple* couple)
{
  currentCouple = couple;
  currentCoupleIndex = couple->GetIndex();
  G4double cross = 0.0;
  if(energy < minPEEnergy) {
    cross = ComputeGeneralLambda(0, 0, idx0, energy);
    peLambda = (thePhotoElectric) 
      ? thePhotoElectric->GetLambda(energy, couple) : 0.0;
    cross += peLambda; 
    
  } else if(energy < minEEEnergy) {
    cross = ComputeGeneralLambda(1, 2, idx1, energy);

  } else if(energy < minMMEnergy) {
    cross = ComputeGeneralLambda(2, 6, idx2, energy);

  } else {
    cross = ComputeGeneralLambda(3, 10, idx3, energy);
  }
  /*
  G4cout << "xs= " << cross << " idxE= " << idxEnergy 
	 << " idxC= " << currentCoupleIndex 
	 << " E= " << energy << G4endl;
  */
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4GammaGeneralProcess::PostStepDoIt(const G4Track& track,
                                                     const G4Step& step)
{
  // In all cases clear number of interaction lengths
  theNumberOfInteractionLengthLeft = -1.0;
  G4double q = G4UniformRand(); 
  G4double x = preStepLambda*G4UniformRand(); 
  G4double p;
  /*
  G4cout << "PostStep: preStepLambda= " << preStepLambda << " PE= " << peLambda
  	 << " q= " << q << G4endl;
  */
  switch (idxEnergy) {
  case 0:
    if(x <= peLambda) { 
      return SampleSecondaries(track, step, thePhotoElectric);
    } else {
      p = GetProbability(1, idx0);
      if(x <= peLambda + (preStepLambda - peLambda)*p) {
	return SampleSecondaries(track, step, theCompton);
      } else if(theRayleigh) {
	return SampleSecondaries(track, step, theRayleigh);
      }
    }
    break;

  case 1:  
    p = GetProbability(3, idx1);
    if(q <= p) {
      return SampleSecondaries(track, step, thePhotoElectric);
    }
    p = GetProbability(4, idx1);
    if(q <= p) {
      return SampleSecondaries(track, step, theCompton);
    }
    p = GetProbability(5, idx1);
    if(q <= p) {
      if(theRayleigh) {
	return SampleSecondaries(track, step, theRayleigh);
      }
    } else if(theGammaNuclear) {
      return SampleSecondaries(track, step, theGammaNuclear);
    }
    break;

  case 2:  
    p = GetProbability(7, idx2);
    if(q <= p) {
      return SampleSecondaries(track, step, theConversionEE);
    }
    p = GetProbability(8, idx2);
    if(q <= p) {
      return SampleSecondaries(track, step, theCompton);
    }
    p = GetProbability(9, idx2);
    if(q <= p) {
      return SampleSecondaries(track, step, thePhotoElectric);
    } else if(theGammaNuclear) {
      return SampleSecondaries(track, step, theGammaNuclear);
    }
    break;

  case 3:  
    p = 1.0 - GetProbability(11, idx3);
    if(q <= p) {
      return SampleSecondaries(track, step, theConversionEE);
    }
    p = 1.0 - GetProbability(12, idx3);
    if(q <= p) {
      return SampleSecondaries(track, step, theCompton);
    }
    p = 1.0 - GetProbability(13, idx3);
    if(q <= p) {
      return SampleSecondaries(track, step, thePhotoElectric);
    } 
    p = 1.0 - GetProbability(14, idx3);
    if(q <= p) {
      if(theGammaNuclear) {
	return SampleSecondaries(track, step, theGammaNuclear);
      }
    } else if(theConversionMM) {
      SelectedProcess(track, theConversionMM);
      return theConversionMM->PostStepDoIt(track, step);
    }
    break;
  }
  // no interaction
  fParticleChange.InitializeForPostStep(track);
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4GammaGeneralProcess::SampleSecondaries(
           const G4Track& track, const G4Step& step, G4HadronicProcess* proc)
{
  SelectedProcess(track, proc);
  proc->GetCrossSectionDataStore()->ComputeCrossSection(track.GetDynamicParticle(),
                                                        track.GetMaterial());
  return theGammaNuclear->PostStepDoIt(track, step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4GammaGeneralProcess::StorePhysicsTable(const G4ParticleDefinition* part,
                                              const G4String& directory,
                                              G4bool ascii)
{
  G4bool yes = true;
  if(!isTheMaster) { return yes; }
  if(thePhotoElectric &&
     !thePhotoElectric->StorePhysicsTable(part, directory, ascii)) 
    { yes = false; }
  if(theCompton && 
     !theCompton->StorePhysicsTable(part, directory, ascii))
    { yes = false; }
  if(theConversionEE && 
     !theConversionEE->StorePhysicsTable(part, directory, ascii))
    { yes = false; }
  if(theRayleigh &&
     !theRayleigh->StorePhysicsTable(part, directory, ascii))
    { yes = false; }

  for(size_t i=0; i<nTables; ++i) {
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
  if(thePhotoElectric &&
     !thePhotoElectric->RetrievePhysicsTable(part, directory, ascii)) 
    { yes = false; }
  if(theCompton && 
     !theCompton->RetrievePhysicsTable(part, directory, ascii))
    { yes = false; }
  if(theConversionEE && 
     !theConversionEE->RetrievePhysicsTable(part, directory, ascii))
    { yes = false; }
  if(theRayleigh &&
     !theRayleigh->RetrievePhysicsTable(part, directory, ascii))
    { yes = false; }

  for(size_t i=0; i<nTables; ++i) {
    if(theT[i]) {
      G4String nam = (0==i || 2==i || 6==i || 10==i) 
	? "LambdaGeneral" + nameT[i] : "ProbGeneral" + nameT[i];
      G4String fnam =  GetPhysicsTableFileName(part,directory,nam,ascii);
      if(!theHandler->RetrievePhysicsTable(i, part, fnam, ascii)) 
	{ yes = false; }
    }
  }
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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
  if(thePhotoElectric) { thePhotoElectric->ProcessDescription(out); }
  if(theCompton)       { theCompton->ProcessDescription(out); }
  if(theConversionEE)  { theConversionEE->ProcessDescription(out); }
  if(theRayleigh)      { theRayleigh->ProcessDescription(out); }
  if(theGammaNuclear)  { theGammaNuclear->ProcessDescription(out); }
  if(theConversionMM)  { theConversionMM->ProcessDescription(out); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4String& G4GammaGeneralProcess::GetProcessName() const
{
  //G4cout << "GetProcessName(): " << selectedProc << G4endl;
  return (selectedProc) ? selectedProc->GetProcessName() 
    : G4VProcess::GetProcessName();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4GammaGeneralProcess::GetProcessSubType() const
{
  return (selectedProc) ? selectedProc->GetProcessSubType() 
    : fGammaGeneralProcess; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
