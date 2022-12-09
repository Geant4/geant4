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
// GEANT4 Class class file
//
//
// File name:     G4VAtomDeexcitation
//
// Author:        Alfonso Mantero & Vladimir Ivanchenko
//
// Creation date: 21.04.2010
//
// Modifications:
//
// Class Description:
//
// Abstract interface to energy loss models

// -------------------------------------------------------------------
//

#include "G4VAtomDeexcitation.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmParameters.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Step.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "Randomize.hh"
#include "G4VParticleChange.hh"
#include "G4EmSecondaryParticleType.hh"
#include "G4Gamma.hh"
#include "G4Log.hh"

#ifdef G4MULTITHREADED
  G4Mutex G4VAtomDeexcitation::atomDeexcitationMutex = G4MUTEX_INITIALIZER;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VAtomDeexcitation::G4VAtomDeexcitation(const G4String& modname) 
  : name(modname)
{
  vdyn.reserve(5);
  theCoupleTable = nullptr;
  gamma = G4Gamma::Gamma();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VAtomDeexcitation::~G4VAtomDeexcitation() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VAtomDeexcitation::InitialiseAtomicDeexcitation()
{
  G4EmParameters* theParameters = G4EmParameters::Instance();
  theParameters->DefineRegParamForDeex(this);

  // Define list of couples
  theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  nCouples = (G4int)theCoupleTable->GetTableSize();

  // needed for unit tests
  std::size_t nn = std::max(nCouples, 1);
  if(activeDeexcitationMedia.size() != nn) {
    activeDeexcitationMedia.resize(nn, false);
    activeAugerMedia.resize(nn, false);
    activePIXEMedia.resize(nn, false);
  }
  if(activeZ.size() != 93) { activeZ.resize(93, false); }

  // initialisation of flags and options
  // normally there is no locksed flags
  if(!isActiveLocked) { isActive  = theParameters->Fluo(); }
  if(!isAugerLocked)  { flagAuger = theParameters->Auger(); }
  if(!isPIXELocked)   { flagPIXE  = theParameters->Pixe(); }
  ignoreCuts = theParameters->DeexcitationIgnoreCut();

  // Define list of regions
  std::size_t nRegions = deRegions.size();
  // check if deexcitation is active for the given run
  if(!isActive && 0 == nRegions) { return; }

  // if no active regions add a world
  if(0 == nRegions) {
    SetDeexcitationActiveRegion("World",isActive,flagAuger,flagPIXE);
    nRegions = deRegions.size();
  }

  if(0 < verbose) {
    G4cout << G4endl;
    G4cout << "### ===  Deexcitation model " << name 
           << " is activated for " << nRegions;
    if(1 == nRegions) { G4cout << " region:" << G4endl; }
    else              { G4cout << " regions:" << G4endl;}
  }

  // Identify active media
  const G4RegionStore* regionStore = G4RegionStore::GetInstance();
  for(std::size_t j=0; j<nRegions; ++j) {
    const G4Region* reg = regionStore->GetRegion(activeRegions[j], false);
    if(nullptr != reg && 0 < nCouples) {
      const G4ProductionCuts* rpcuts = reg->GetProductionCuts();
      if(0 < verbose) {
        G4cout << "          " << activeRegions[j]
               << "  " << deRegions[j]  << "  " << AugerRegions[j]
               << "  " << PIXERegions[j] << G4endl;  
      }
      for(G4int i=0; i<nCouples; ++i) {
        const G4MaterialCutsCouple* couple =
          theCoupleTable->GetMaterialCutsCouple(i);
        if (couple->GetProductionCuts() == rpcuts) {
          activeDeexcitationMedia[i] = deRegions[j];
          activeAugerMedia[i] = AugerRegions[j];
          activePIXEMedia[i] = PIXERegions[j];
        }
      }
    }
  }
  std::size_t nelm = G4Element::GetNumberOfElements();
  //G4cout << nelm << G4endl;
  for(std::size_t k=0; k<nelm; ++k) {
    G4int Z = (*(G4Element::GetElementTable()))[k]->GetZasInt();
    if(Z > 5 && Z < 93) { 
      activeZ[Z] = true;
      //G4cout << "!!! Active de-excitation Z= " << Z << G4endl;  
    }
  }

  // Initialise derived class
  InitialiseForNewRun();

  if(0 < verbose && flagAuger) {
    G4cout << "### ===  Auger flag: " << flagAuger 
           << G4endl;
  }
  if(0 < verbose) {
    G4cout << "### ===  Ignore cuts flag:   " << ignoreCuts
           << G4endl;
  }
  if(0 < verbose && flagPIXE) {
    G4cout << "### ===  PIXE model for hadrons: " 
           << theParameters->PIXECrossSectionModel()
           << G4endl;  
    G4cout << "### ===  PIXE model for e+-:     " 
           << theParameters->PIXEElectronCrossSectionModel()
           << G4endl;  
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4VAtomDeexcitation::SetDeexcitationActiveRegion(const G4String& rname,
                                                 G4bool valDeexcitation,
                                                 G4bool valAuger,
                                                 G4bool valPIXE)
{
  // no PIXE in parallel world
  if(rname == "DefaultRegionForParallelWorld") { return; }

  G4String ss = rname;
  /*  
  G4cout << "### G4VAtomDeexcitation::SetDeexcitationActiveRegion " << ss 
         << "  " << valDeexcitation << "  " << valAuger
         << "  " << valPIXE << G4endl;
  */
  if(ss == "world" || ss == "World" || ss == "WORLD") {
    ss = "DefaultRegionForTheWorld";
  }
  std::size_t n = deRegions.size();
  for(std::size_t i=0; i<n; ++i) {
 
    // Region already exist
    if(ss == activeRegions[i]) {
      deRegions[i] = valDeexcitation;
      AugerRegions[i] = valAuger;
      PIXERegions[i] = valPIXE;
      return;  
    }
  }
  // New region
  activeRegions.push_back(ss);
  deRegions.push_back(valDeexcitation);
  AugerRegions.push_back(valAuger);
  PIXERegions.push_back(valPIXE);

  // if de-excitation defined for the world volume 
  // it should be active for all G4Regions
  if(ss == "DefaultRegionForTheWorld") {
    G4RegionStore* regions = G4RegionStore::GetInstance();
    std::size_t nn = regions->size();
    for(std::size_t i=0; i<nn; ++i) {
      if(ss == (*regions)[i]->GetName()) { continue; }
      SetDeexcitationActiveRegion((*regions)[i]->GetName(), valDeexcitation,
                                  valAuger, valPIXE);
                                  
    }
  }
}

void G4VAtomDeexcitation::GenerateParticles(std::vector<G4DynamicParticle*>* v,
                                            const G4AtomicShell* as, 
                                            G4int Z, G4int idx)
{
  G4double gCut = DBL_MAX;
  if(ignoreCuts) {
    gCut = 0.0;
  } else if (nullptr != theCoupleTable) {
    gCut = (*(theCoupleTable->GetEnergyCutsVector(0)))[idx];
  }
  if(gCut < as->BindingEnergy()) {
    G4double eCut = DBL_MAX;
    if(CheckAugerActiveRegion(idx)) {
      if(ignoreCuts) {
        eCut = 0.0;
      } else if (nullptr != theCoupleTable) {
        eCut = (*(theCoupleTable->GetEnergyCutsVector(1)))[idx];
      }
    }
    GenerateParticles(v, as, Z, gCut, eCut);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4VAtomDeexcitation::AlongStepDeexcitation(std::vector<G4Track*>& tracks,
                                           const G4Step& step, 
                                           G4double& eLossMax,
                                           G4int coupleIndex)
{
  G4double truelength = step.GetStepLength();
  if(!flagPIXE && !activePIXEMedia[coupleIndex]) { return; }
  if(eLossMax <= 0.0 || truelength <= 0.0)       { return; }

  // step parameters
  const G4StepPoint* preStep = step.GetPreStepPoint();
  const G4ThreeVector prePos = preStep->GetPosition();
  const G4ThreeVector delta = step.GetPostStepPoint()->GetPosition() - prePos;
  const G4double preTime = preStep->GetGlobalTime();
  const G4double dt = step.GetPostStepPoint()->GetGlobalTime() - preTime;

  // particle parameters
  const G4Track* track = step.GetTrack();
  const G4ParticleDefinition* part = track->GetDefinition();
  G4double ekin = preStep->GetKineticEnergy(); 

  // media parameters
  G4double gCut = (*theCoupleTable->GetEnergyCutsVector(0))[coupleIndex];
  if(ignoreCuts) { gCut = 0.0; }
  G4double eCut = DBL_MAX;
  if(CheckAugerActiveRegion(coupleIndex)) { 
    eCut = (*theCoupleTable->GetEnergyCutsVector(1))[coupleIndex];
    if(ignoreCuts) { eCut = 0.0; }    
  }

  //G4cout<<"!Sample PIXE gCut(MeV)= "<<gCut<<"  eCut(MeV)= "<<eCut
  //        <<" Ekin(MeV)= " << ekin/MeV << G4endl;

  const G4Material* material = preStep->GetMaterial();
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = 
    material->GetVecNbOfAtomsPerVolume();
  const std::size_t nelm = material->GetNumberOfElements();

  // loop over deexcitations
  for(std::size_t i=0; i<nelm; ++i) {
    G4int Z = (*theElementVector)[i]->GetZasInt();
    if(activeZ[Z] && Z < 93) {  
      G4int nshells = 
        std::min(9,(*theElementVector)[i]->GetNbOfAtomicShells());
      G4double rho = truelength*theAtomNumDensityVector[i];
      //G4cout<<"   Z "<< Z <<" is active  x(mm)= " << truelength/mm << G4endl;
      for(G4int ii=0; ii<nshells; ++ii) {
        auto as = (G4AtomicShellEnumerator)(ii);
        const G4AtomicShell* shell = GetAtomicShell(Z, as);
        const G4double bindingEnergy = shell->BindingEnergy();

        if(gCut > bindingEnergy) { break; }

        if(eLossMax > bindingEnergy) { 
          G4double sig = rho*
            GetShellIonisationCrossSectionPerAtom(part, Z, as, ekin, material);

          // mfp is mean free path in units of step size
          if(sig > 0.0) {
            G4double mfp = 1.0/sig;
            G4double stot = 0.0;
            //G4cout << " Shell " << ii << " mfp(mm)= " << mfp/mm << G4endl;
            // sample ionisation points
            do {
              stot -= mfp*G4Log(G4UniformRand());
              if( stot > 1.0 || eLossMax < bindingEnergy) { break; }
              // sample deexcitation
              vdyn.clear();
              GenerateParticles(&vdyn, shell, Z, gCut, eCut); 
              std::size_t nsec = vdyn.size();
              if(nsec > 0) {
                G4ThreeVector r = prePos  + stot*delta;
                G4double time   = preTime + stot*dt;
                for(std::size_t j=0; j<nsec; ++j) {
                  G4DynamicParticle* dp = vdyn[j];
                  G4double e = dp->GetKineticEnergy();

                  // save new secondary if there is enough energy
                  if(eLossMax >= e) {
                    eLossMax -= e;                    
                    G4Track* t = new G4Track(dp, time, r);

                    // defined secondary type
                    if(dp->GetDefinition() == gamma) { 
                      t->SetCreatorModelID(_GammaPIXE);
                    } else {
                      t->SetCreatorModelID(_ePIXE);
                    }
                    tracks.push_back(t);
                  } else {
                    delete dp;
                  }
                }
              }
              // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
            } while (stot < 1.0);
          }
        }
      }
    } 
  }
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
