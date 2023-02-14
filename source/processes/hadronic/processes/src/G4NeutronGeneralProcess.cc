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
// File name:     G4NeutronGeneralProcess
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 08.08.2022
//
// Modifications:
//
// Class Description:
//

// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4NeutronGeneralProcess.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicProcess.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4VParticleChange.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4HadronicParameters.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4Neutron.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronElasticXS.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4Threading.hh"

#include "G4Log.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4HadDataHandler* G4NeutronGeneralProcess::theHandler = nullptr;

G4String G4NeutronGeneralProcess::nameT[nTables] = {"0","1","2","3","4"};

G4NeutronGeneralProcess::G4NeutronGeneralProcess(const G4String& pname)
: G4HadronicProcess(pname),
  fMinEnergy(1*CLHEP::keV),
  fMiddleEnergy(20*CLHEP::MeV),
  fMaxEnergy(100*CLHEP::TeV),
  fTimeLimit(10*CLHEP::microsecond)
{
  SetVerboseLevel(1);
  SetProcessSubType(fNeutronGeneral);

  fNeutron = G4Neutron::Neutron();

  if(G4Threading::IsWorkerThread()) {
    isMaster = false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4NeutronGeneralProcess::~G4NeutronGeneralProcess()
{
  if(isMaster) {
    delete theHandler;
    theHandler = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4NeutronGeneralProcess::IsApplicable(const G4ParticleDefinition&)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4NeutronGeneralProcess::SetInelasticProcess(G4HadronicProcess* ptr)
{
  fInelastic = ptr;
  fXSSInelastic = ptr->GetCrossSectionDataStore();
  fInelasticXS = InitialisationXS(ptr);
  if(nullptr == fInelasticXS) {
    fInelasticXS = new G4NeutronInelasticXS();
    ptr->AddDataSet(fInelasticXS);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4NeutronGeneralProcess::SetElasticProcess(G4HadronicProcess* ptr)
{
  fElastic = ptr;
  fXSSElastic = ptr->GetCrossSectionDataStore();
  fElasticXS = InitialisationXS(ptr);
  if(nullptr == fElasticXS) {
    fElasticXS = new G4NeutronElasticXS();
    ptr->AddDataSet(fElasticXS);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4NeutronGeneralProcess::SetCaptureProcess(G4HadronicProcess* ptr)
{
  fCapture = ptr;
  fXSSCapture = ptr->GetCrossSectionDataStore();
  fCaptureXS = InitialisationXS(ptr);
  if(nullptr == fCaptureXS) {
    fCaptureXS = new G4NeutronCaptureXS();
    ptr->AddDataSet(fCaptureXS);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VCrossSectionDataSet*
G4NeutronGeneralProcess::InitialisationXS(G4HadronicProcess* proc)
{
  G4VCrossSectionDataSet* ptr = nullptr;
  auto xsv = proc->GetCrossSectionDataStore()->GetDataSetList();
  if(!xsv.empty()) {
    ptr = xsv[0];
  }
  return ptr;
}

//....Ooooo0ooooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4NeutronGeneralProcess::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  if(1 < verboseLevel) {
    G4cout << "G4NeutronGeneralProcess::PreparePhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
 	   << " isMaster: " << isMaster << G4endl;
  }
  G4bool noEl = (nullptr == fElastic);
  G4bool noInel = (nullptr == fInelastic);
  G4bool noCap = (nullptr == fCapture);
  if(noEl || noInel || noCap) {
    G4ExceptionDescription ed;
    ed << "Incomplete configuration of the neutron general process." << G4endl;
    if(noEl) { ed << "Neutron elastic process is not defined" << G4endl; }
    if(noInel) { ed << "Neutron inelastic process is not defined" << G4endl; }
    if(noCap) { ed << "Neutron capture process is not defined" << G4endl; }
    G4Exception ("G4NeutronGeneralProcess::PreparePhysicsTable(..)", "had001",
                 FatalException, ed, "");
    return;
  }

  G4HadronicParameters* param = G4HadronicParameters::Instance();

  SetVerboseLevel(param->GetVerboseLevel());
  fMaxEnergy = std::max(100*MeV, param->GetMaxEnergy());
  if(param->ApplyFactorXS()) {
    fXSFactorEl = param->XSFactorNucleonElastic();
    fXSFactorInel = param->XSFactorNucleonInelastic();
  }

  fElastic->PreparePhysicsTable(part);
  fInelastic->PreparePhysicsTable(part);
  fCapture->PreparePhysicsTable(part);

  std::size_t nmat = G4Material::GetNumberOfMaterials();
  G4MaterialTable* matTable = G4Material::GetMaterialTable();

  std::size_t nmax = 0;
  for(std::size_t i=0; i<nmat; ++i) {
    std::size_t nelm = (*matTable)[i]->GetNumberOfElements();
    nmax = std::max(nmax, nelm);
  }
  fXsec.resize(nmax);

  if(isMaster) {
    if(nullptr == theHandler) {
      theHandler = new G4HadDataHandler(nTables);
    }

    fMaxEnergy = std::max(fMaxEnergy, param->GetMaxEnergy());
    nLowE *= G4lrint(std::log10(fMiddleEnergy/fMinEnergy));
    nHighE *= G4lrint(std::log10(fMaxEnergy/fMiddleEnergy));

    G4PhysicsVector* vec = nullptr;
    G4PhysicsLogVector aVector(fMinEnergy, fMiddleEnergy, nLowE, false);
    G4PhysicsLogVector bVector(fMiddleEnergy, fMaxEnergy, nHighE, false);

    for(std::size_t i=0; i<nTables; ++i) {
      G4PhysicsTable* table = new G4PhysicsTable();
      theHandler->UpdateTable(table, i);
      table->resize(nmat);
      for(std::size_t j=0; j<nmat; ++j) {
	vec = (*table)[j];
	if (nullptr == vec) {
	  if(i <= 2) {
	    vec = new G4PhysicsVector(aVector);
	  } else {
	    vec = new G4PhysicsVector(bVector);
	  }
	  G4PhysicsTableHelper::SetPhysicsVector(table, j, vec);
	}
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4NeutronGeneralProcess::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  if(1 < verboseLevel) {
    G4cout << "### G4NeutronGeneralProcess::BuildPhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << G4endl;
  }
  fElastic->BuildPhysicsTable(part);
  fInelastic->BuildPhysicsTable(part);
  fCapture->BuildPhysicsTable(part);

  if(isMaster) {
    std::size_t nmat = G4Material::GetNumberOfMaterials();
    G4MaterialTable* matTable = G4Material::GetMaterialTable();

    auto tables = theHandler->GetTables();

    G4double sigEl(0.), sigInel(0.), sigCap(0.), val(0.), sum(0.);

    for(std::size_t i=0; i<nmat; ++i) {
      const G4Material* mat = (*matTable)[i];

      // energy interval 0
      std::size_t nn = (*(tables[0]))[i]->GetVectorLength();
      if(1 < verboseLevel) {
	G4cout << "======= Zone 0 ======= N= " << nn
	       << " for " << mat->GetName() << G4endl;
      }
      for(std::size_t j=0; j<nn; ++j) {
	G4double e = (*(tables[0]))[i]->Energy(j);
	G4double loge = G4Log(e);
        sigEl = fXSFactorEl*ComputeCrossSection(fElasticXS, mat, e, loge);
        sigInel = fXSFactorInel*ComputeCrossSection(fInelasticXS, mat, e, loge);
        sigCap = ComputeCrossSection(fCaptureXS, mat, e, loge);
	sum = sigEl + sigInel + sigCap;
	if(1 < verboseLevel) {
	  G4cout << j << ". E= " << e << " xs=" << sum << " sigEl=" << sigEl
		 << " sigInel=" << sigInel << " sigCap=" << sigCap << G4endl;
	}
	(*(tables[0]))[i]->PutValue(j, sum);
	val = sigEl/sum;
	(*(tables[1]))[i]->PutValue(j, val);
	val = (sigEl + sigInel)/sum;
	(*(tables[2]))[i]->PutValue(j, val);
      }

      // energy interval 1
      nn = (*(tables[3]))[0]->GetVectorLength();
      if(1 < verboseLevel) {
	G4cout << "======= Zone 1 ======= N= " << nn << G4endl;
      }
      for(std::size_t j=0; j<nn; ++j) {
	G4double e = (*(tables[3]))[i]->Energy(j);
	G4double loge = G4Log(e);
	sigEl = fXSFactorEl*ComputeCrossSection(fElasticXS, mat, e, loge);
	sigInel = fXSFactorInel*ComputeCrossSection(fInelasticXS, mat, e, loge);
	sum = sigEl + sigInel;
	if(1 < verboseLevel) {
	  G4cout << j << ". E= " << e << " xs=" << sum << " sigEl=" << sigEl
		 << " sigInel=" << sigInel << " factInel=" << fXSFactorInel 
                 << G4endl;
	}
	(*(tables[3]))[i]->PutValue(j, sum);
	val = sigInel/sum;
	(*(tables[4]))[i]->PutValue(j, val);
      }
    }
  }
  if(1 < verboseLevel) {
    G4cout << "### G4VEmProcess::BuildPhysicsTable() done for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4NeutronGeneralProcess::ComputeCrossSection(G4VCrossSectionDataSet* xs,
                                             const G4Material* mat,
                                             G4double e, G4double loge)
{
  const G4double* natom = mat->GetVecNbOfAtomsPerVolume();
  G4int nelm = (G4int)mat->GetNumberOfElements();
  G4double sig = 0.0;
  for(G4int i=0; i<nelm; ++i) {
    sig += natom[i]*xs->ComputeCrossSectionPerElement(e, loge, fNeutron,
                                                      mat->GetElement(i), mat);
  }
  return sig;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4NeutronGeneralProcess::StartTracking(G4Track*)
{
  theNumberOfInteractionLengthLeft = -1.0;
  fCurrMat = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4NeutronGeneralProcess::PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double previousStepSize,
                             G4ForceCondition* condition)
{
  *condition = NotForced;

  // time limit
  if(track.GetGlobalTime() >= fTimeLimit) {
    fLambda = 0.0;
    return 0.0;
  }

  // recompute total cross section if needed
  CurrentCrossSection(track);

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

  G4double x = theNumberOfInteractionLengthLeft * currentInteractionLength;
  /*
  G4cout << "PostStepGetPhysicalInteractionLength: e= " << energy
	 << " idxe= " << idxEnergy << "  xs= " << fLambda
	 << " x= " << x << G4endl;
  */
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4NeutronGeneralProcess::PostStepDoIt(const G4Track& track,
                                                         const G4Step& step)
{
  fSelectedProc = this;
  // time limit
  if(0.0 == fLambda) {
    theTotalResult->Initialize(track);
    theTotalResult->ProposeTrackStatus(fStopAndKill);
    return theTotalResult;
  }
  // In all cases clear number of interaction lengths
  theNumberOfInteractionLengthLeft = -1.0;
  G4double q = G4UniformRand();
  /*  
  G4cout << "PostStep: preStepLambda= " << fLambda << " idxE= " << idxEnergy
         << " matIndex=" << matIndex << G4endl;
  */
  if (0 == idxEnergy) {
    if(q <= GetProbability(1)) {
      SelectedProcess(step, fElastic, fXSSElastic);
    } else if(q <= GetProbability(2)) {
      SelectedProcess(step, fInelastic, fXSSInelastic);
    } else {
      SelectedProcess(step, fCapture, fXSSCapture);
    }
  } else {
    if(q <= GetProbability(4)) {
      SelectedProcess(step, fInelastic, fXSSInelastic);
    } else {
      SelectedProcess(step, fElastic, fXSSElastic);
    }
  }
  // total cross section is needed for selection of an element
  if(fCurrMat->GetNumberOfElements() > 1) {
    fCurrentXSS->ComputeCrossSection(track.GetDynamicParticle(), fCurrMat);
  }
  /*
    G4cout << "## neutron E(MeV)=" << fCurrE << " inside " << fCurrMat->GetName() 
	 << fSelectedProc->GetProcessName()
	 << " time(ns)=" << track.GetGlobalTime()/ns << G4endl; 
  */
  // sample secondaries
  return fSelectedProc->PostStepDoIt(track, step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool
G4NeutronGeneralProcess::StorePhysicsTable(const G4ParticleDefinition* part,
                                           const G4String& directory,
                                           G4bool ascii)
{
  G4bool yes = true;
  if(!isMaster) { return yes; }
  for(std::size_t i=0; i<nTables; ++i) {
    G4String nam = (0==i || 3==i)
	? "LambdaNeutronGeneral" + nameT[i] : "ProbNeutronGeneral" + nameT[i];
    G4String fnam =  GetPhysicsTableFileName(part, directory, nam, ascii);
    auto table = theHandler->Table(i);
    if(nullptr == table || !table->StorePhysicsTable(fnam, ascii)) { 
      yes = false;
    }
  }
  return yes;
}

//....Ooooo0ooooo ........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4NeutronGeneralProcess::GetMeanFreePath(const G4Track& track,
						  G4double,
						  G4ForceCondition* condition)
{
  *condition = NotForced;
  // recompute total cross section if needed
  CurrentCrossSection(track);
  return currentInteractionLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4NeutronGeneralProcess::ProcessDescription(std::ostream& out) const
{
  fElastic->ProcessDescription(out);
  fInelastic->ProcessDescription(out);
  fCapture->ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4String& G4NeutronGeneralProcess::GetSubProcessName() const
{
  return (nullptr != fSelectedProc) ? fSelectedProc->GetProcessName()
    : G4VProcess::GetProcessName();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4NeutronGeneralProcess::GetSubProcessSubType() const
{
  return (nullptr != fSelectedProc) ? fSelectedProc->GetProcessSubType()
    : fNeutronGeneral;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4VProcess* G4NeutronGeneralProcess::GetCreatorProcess() const
{
  return fSelectedProc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
