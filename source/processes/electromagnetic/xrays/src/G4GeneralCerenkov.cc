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
// G4GeneralCerenkov
//
// Created 25.05.2025 V.Ivanchenko
//
// --------------------------------------------------------------------

#include "G4GeneralCerenkov.hh"
#include "G4StandardCerenkovModel.hh"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalParameters.hh"
#include "G4PhysicsModelCatalog.hh"
#include "G4EmProcessSubType.hh"
#include "G4LogicalVolumeStore.hh"

std::vector<std::vector<const G4LogicalVolume*>* >* G4GeneralCerenkov::fLV = nullptr;
std::vector<G4VXRayModel*>* G4GeneralCerenkov::fSharedModels = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4GeneralCerenkov::G4GeneralCerenkov(const G4String& nam, G4ProcessType type)
  : G4VDiscreteProcess(nam, type)
{
  secID = G4PhysicsModelCatalog::GetModelID("model_Cerenkov");
  SetProcessSubType(fCerenkov);
  if (nullptr == fLV) {
    // initialise static data    
    fSharedModels = new std::vector<G4VXRayModel*>;
    fLV = new std::vector<std::vector<const G4LogicalVolume*>* >;
    fLVNames = new std::vector<G4String>;

    isInitializer = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4GeneralCerenkov::~G4GeneralCerenkov()
{
  if (isInitializer) {
    delete fSharedModels;
    fSharedModels = nullptr;
    delete fLVNames;
    for (auto const & p : *fLV) {
      delete p;
    }
    delete fLV;
    fLV = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool G4GeneralCerenkov::IsApplicable(const G4ParticleDefinition&)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4GeneralCerenkov::AddModelForVolume(G4VXRayModel* model,
					  const G4String& nameLogVolume)
{
  if (isPrepared || !isInitializer || nullptr == model) {
     G4ExceptionDescription ed;
     G4String nam;
     if (model != nullptr) { nam = model->GetName(); }
     ed << " Attempt to add Cerenkov model <" << nam << "> for LogicalVolume "
	<< nameLogVolume << " is failed!\n isPrepared:" << isPrepared
	<< " isInitilizer:" << isInitializer;
     G4Exception("G4GeneralCerenkov::AddModelForVolume", "em0304",
		 FatalException, ed, "");
     return;
  }
  fSharedModels->push_back(model);
  fLVNames->push_back(nameLogVolume);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4GeneralCerenkov::PreparePhysicsTable(const G4ParticleDefinition&)
{
  // definition of models is done only once
  if (isPrepared) { return; }
  isPrepared = true;
  
  const G4OpticalParameters* params = G4OpticalParameters::Instance();
  fMaxBetaChange = params->GetCerenkovMaxBetaChange();
  fMaxPhotons = params->GetCerenkovMaxPhotonsPerStep();
  fStackingFlag = params->GetCerenkovStackPhotons();
  fTrackSecondariesFirst = params->GetCerenkovTrackSecondariesFirst();
  verboseLevel = params->GetCerenkovVerboseLevel();

  auto nmod = fSharedModels->size();
  if (0 == nmod) {
    // the default model is added without association with a logical volume
    G4VXRayModel* mod = new G4StandardCerenkovModel();    
    fSharedModels->push_back(mod);
    nmod = 1;
  }

  fSecondaries.reserve(fMaxPhotons);
  nModels = (G4int)nmod;

  // fill static data structures
  if (isInitializer) {
    const G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
    const auto & modAndVol = params->ActiveVolumes();

    // preparation of logical volume vector per model
    fLV->reserve(nmod);
    for (G4int i=0; i<nModels; ++i) {
      auto v = new std::vector<const G4LogicalVolume*>;
      fLV->push_back(v);
    }
    
    for (auto const & lv : *lvs) {
      // only volumes with material property defined are considered
      auto const MPT = lv->GetMaterial()->GetMaterialPropertiesTable();
      if (nullptr == MPT) { continue; }
      
      const G4String& lvname = lv->GetName();
      G4bool ok{false};
      // search for the default model in the list
      if (!modAndVol.empty()) {
	for (auto const & it : modAndVol) {
	  if (it.second == lvname) {
	    if (kCerenkovDefault == it.first) {
	      (*fLV)[0]->push_back(lv);
	      fLVNames->push_back(lvname);
	      ok = true;
	      break;
	    }
	  }
	}
      }
      if (ok) { continue; }

      // search in external models
      for (G4int i=0; i<nModels; ++i) {
        if (lvname == (*fLVNames)[i]) {
	  (*fLV)[i]->push_back(lv);
	  ok = true;
	  break;
	}
      }
      if (ok) { continue; }

      // temporary for backward compatibility search for a RINDEX of the volume
      if (nullptr != MPT->GetProperty(kRINDEX)) {
	(*fLV)[0]->push_back(lv);
	fLVNames->push_back(lvname);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4GeneralCerenkov::BuildPhysicsTable(const G4ParticleDefinition&)
{
  // Initialisation of models is done only once
  if (isBuilt) { return; }
  isBuilt = true;
  if (nModels == 0) { return; }

  // worker thread
  if (!isInitializer) {
    // worker initialisation - clone master models
    fModels.reserve(nModels);
    for (G4int i=0; i<nModels; ++i) {
      auto newmod = new G4VXRayModel(*((*fSharedModels)[i]));
      fModels.push_back(newmod);
      G4double b = newmod->Initialise((*fLV)[i]);
      fBetaMin = std::min(fBetaMin, b);
    }
  } else if (verboseLevel > 0) {
    // needed for printout
    std::size_t nn = 0;
    for (G4int i=0; i<nModels; ++i) {
      G4double b = (*fSharedModels)[i]->Initialise((*fLV)[i]);
      fBetaMin = std::min(fBetaMin, b);
      nn += ((*fLV)[i])->size();
    }
    G4long pres = G4cout.precision();
    G4cout.precision(6);
    G4cout << "  " << GetProcessName() << std::setw(20) << "  fBetaMin=" << fBetaMin
	   << "  fMaxBetaChange=" << fMaxBetaChange << G4endl;
    G4cout << std::setw(20) << "fMaxNphot=" << fMaxPhotons
	   << "  Nlv=" << nn << "  fStackingFlag:" << fStackingFlag
	   << "  fTrackSecondariesFirst:" << fTrackSecondariesFirst
	   << G4endl;
    for (G4int i=0; i<nModels; ++i) {
      G4int n = (G4int)((*fLV)[i]->size());
      G4cout << std::setw(10) << (*fSharedModels)[i]->GetName()
	     << std::setw(30) << "Nvolumes=" << n << "  Volumes:" << G4endl;
      G4cout << std::setw(12);
      for (G4int j=0; j<n; ++j) {
	G4cout << (*((*fLV)[i]))[j]->GetName() << " ";
	if (0 != j && (j/5)*5 == j) {
	  G4cout << G4endl << std::setw(12);
	}
      }
      G4cout << G4endl;
    }
    G4cout.precision(pres);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double
G4GeneralCerenkov::PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
							G4double,
							G4ForceCondition* cond)
{
  *cond = NotForced;
  G4double limit = DBL_MAX;
  auto const dp = aTrack.GetDynamicParticle();
  if (dp->GetCharge() == 0.0) { return limit; }
  fCurrentModel = nullptr;
  fPreStepBeta = dp->GetBeta();
  if (fPreStepBeta <= fBetaMin) { return limit; }

  auto volume = aTrack.GetVolume();
  if (nullptr == volume) { return limit; }

  fCurrentLV = volume->GetLogicalVolume();
  auto const MPT = fCurrentLV->GetMaterial()->GetMaterialPropertiesTable();
  if (nullptr == MPT) { return limit; }

  G4bool ok{false};
  for (G4int i=0; i<nModels; ++i) {
    auto const v = (*fLV)[i];
    std::size_t nn = v->size();
    for (std::size_t j = 0; j < nn; ++j) {
      if ((*v)[j] == fCurrentLV) {
	fCurrentModel = fModels[i];
	if (fCurrentModel->StepLimit(j, aTrack, fPreStepBeta, limit)) {
	  *cond = StronglyForced;
	}
	ok = true;
	break;
      }
    }
    if (ok) { break; }
  }
  return limit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VParticleChange* G4GeneralCerenkov::PostStepDoIt(const G4Track& aTrack,
						   const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);
  if (fCurrentModel == nullptr) { return &aParticleChange; }
  fCurrentModel->SampleXRays(fSecondaries, aStep);
  if (!fSecondaries.empty()) {

    // X-rays
    auto touch = aStep.GetPreStepPoint()->GetTouchableHandle();
    G4int parent = aTrack.GetTrackID();
    for (auto & t : fSecondaries) {
      t->SetTouchableHandle(touch);
      t->SetParentID(parent);
      t->SetCreatorModelID(secID);
      aParticleChange.AddSecondary(t);
    }
    fSecondaries.clear();

    // primary track suspended
    if (fTrackSecondariesFirst && aTrack.GetTrackStatus() == fAlive) {
      aParticleChange.ProposeTrackStatus(fSuspend);
    }
  }
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4GeneralCerenkov::ProcessDescription(std::ostream& out) const
{
  out << "The Cerenkov effect simulates optical photons created by the\n";
  out << "passage of charged particles through matter. Materials need\n";
  out << "to have the property RINDEX (refractive index) defined." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4GeneralCerenkov::SetTrackSecondariesFirst(const G4bool state)
{
  fTrackSecondariesFirst = state;
  G4OpticalParameters::Instance()->SetCerenkovTrackSecondariesFirst(
    fTrackSecondariesFirst);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4GeneralCerenkov::SetMaxBetaChangePerStep(const G4double value)
{
  fMaxBetaChange = value;
  G4OpticalParameters::Instance()->SetCerenkovMaxBetaChange(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4GeneralCerenkov::SetMaxNumPhotonsPerStep(const G4int NumPhotons)
{
  fMaxPhotons = NumPhotons;
  G4OpticalParameters::Instance()->SetCerenkovMaxPhotonsPerStep(fMaxPhotons);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4GeneralCerenkov::SetStackPhotons(const G4bool stackingFlag)
{
  fStackingFlag = stackingFlag;
  G4OpticalParameters::Instance()->SetCerenkovStackPhotons(fStackingFlag);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4GeneralCerenkov::SetVerboseLevel(G4int verbose)
{
  verboseLevel = verbose;
  G4OpticalParameters::Instance()->SetCerenkovVerboseLevel(verboseLevel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double
G4GeneralCerenkov::GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*)
{
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
