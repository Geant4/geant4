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
// $Id$
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VAtomDeexcitation::G4VAtomDeexcitation(const G4String& modname, 
					 const G4String& pname) 
  : lowestKinEnergy(keV), verbose(1), name(modname), namePIXE(pname), 
    nameElectronPIXE(""), isActive(false), flagAuger(false), flagPIXE(false)
{
  vdyn.reserve(5);
  theCoupleTable = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VAtomDeexcitation::~G4VAtomDeexcitation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VAtomDeexcitation::InitialiseAtomicDeexcitation()
{
  // Define list of couples
  theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  // needed for unit tests
  if(0 == numOfCouples) { numOfCouples = 1; }

  activeDeexcitationMedia.resize(numOfCouples, false);
  activeAugerMedia.resize(numOfCouples, false);
  activePIXEMedia.resize(numOfCouples, false);
  activeZ.resize(93, false);

  // check if deexcitation is active for the given run
  if( !isActive ) { return; }

  // Define list of regions
  size_t nRegions = deRegions.size();

  if(0 == nRegions) {
    SetDeexcitationActiveRegion("World",isActive,flagAuger,flagPIXE);
    nRegions = 1;
  }

  if(0 < verbose) {
    G4cout << G4endl;
    G4cout << "### ===  Deexcitation model " << name 
	   << " is activated for " << nRegions;
    if(1 == nRegions) { G4cout << " region:" << G4endl; }
    else              { G4cout << " regions:" << G4endl;}
  }

  // Identify active media
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  for(size_t j=0; j<nRegions; ++j) {
    const G4Region* reg = regionStore->GetRegion(activeRegions[j], false);
    const G4ProductionCuts* rpcuts = reg->GetProductionCuts();
    if(0 < verbose) {
      G4cout << "          " << activeRegions[j] << G4endl;  
    }
  
    for(size_t i=0; i<numOfCouples; ++i) {
      const G4MaterialCutsCouple* couple =
	theCoupleTable->GetMaterialCutsCouple(i);
      if (couple->GetProductionCuts() == rpcuts) {
	activeDeexcitationMedia[i] = deRegions[j];
	activeAugerMedia[i] = AugerRegions[j];
	activePIXEMedia[i] = PIXERegions[j];
	const G4Material* mat = couple->GetMaterial();
	const G4ElementVector* theElementVector = 
	  mat->GetElementVector();
	G4int nelm = mat->GetNumberOfElements();
	if(deRegions[j]) {
	  for(G4int k=0; k<nelm; ++k) {
	    G4int Z = G4lrint(((*theElementVector)[k])->GetZ());
	    if(Z > 5 && Z < 93) { 
	      activeZ[Z] = true;
	      //G4cout << "!!! Active de-excitation Z= " << Z << G4endl;  
	    }
	  }
	}
      }
    }
  }

  // Initialise derived class
  InitialiseForNewRun();

  if(0 < verbose && flagPIXE) {
    G4cout << "### ===  PIXE model for hadrons: " << namePIXE
	   << "  " <<  IsPIXEActive()
	   << G4endl;  
    G4cout << "### ===  PIXE model for e+-:     " << nameElectronPIXE
	   << "  " <<  IsPIXEActive()
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
  G4String ss = rname;
  //G4cout << "### G4VAtomDeexcitation::SetDeexcitationActiveRegion " << ss 
  //	 << "  " << valDeexcitation << "  " << valAuger
  //	 << "  " << valPIXE << G4endl;
  if(ss == "world" || ss == "World" || ss == "WORLD") {
    ss = "DefaultRegionForTheWorld";
  }
  size_t n = deRegions.size();
  if(n > 0) {
    for(size_t i=0; i<n; ++i) {
 
      // Region already exist
      if(ss == activeRegions[i]) {
	deRegions[i] = valDeexcitation;
	AugerRegions[i] = valAuger;
	PIXERegions[i] = valPIXE;
	return; 
      } 
    }
  }
  // New region
  activeRegions.push_back(ss);
  deRegions.push_back(valDeexcitation);
  AugerRegions.push_back(valAuger);
  PIXERegions.push_back(valPIXE);
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
  G4ThreeVector prePos = preStep->GetPosition();
  G4ThreeVector delta = step.GetPostStepPoint()->GetPosition() - prePos;
  G4double preTime = preStep->GetGlobalTime();
  G4double dt = step.GetPostStepPoint()->GetGlobalTime() - preTime;

  // particle parameters
  const G4Track* track = step.GetTrack();
  const G4ParticleDefinition* part = track->GetDefinition();
  G4double ekin = preStep->GetKineticEnergy(); 

  // media parameters
  G4double gCut = (*theCoupleTable->GetEnergyCutsVector(0))[coupleIndex];
  G4double eCut = DBL_MAX;
  if(CheckAugerActiveRegion(coupleIndex)) { 
    eCut = (*theCoupleTable->GetEnergyCutsVector(1))[coupleIndex];
  }

  //G4cout<<"!Sample PIXE gCut(MeV)= "<<gCut<<"  eCut(MeV)= "<<eCut
  //	<<" Ekin(MeV)= " << ekin/MeV << G4endl;

  const G4Material* material = preStep->GetMaterial();
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
  G4int nelm = material->GetNumberOfElements();

  // loop over deexcitations
  for(G4int i=0; i<nelm; ++i) {
    G4int Z = G4lrint((*theElementVector)[i]->GetZ());
    if(activeZ[Z] && Z < 93) {  
      G4int nshells = std::min(9,(*theElementVector)[i]->GetNbOfAtomicShells());
      G4double rho = truelength*theAtomNumDensityVector[i];
      //G4cout << "   Z " << Z <<" is active  x(mm)= " << truelength/mm << G4endl;
      for(G4int ii=0; ii<nshells; ++ii) {
	G4AtomicShellEnumerator as = G4AtomicShellEnumerator(ii);
	const G4AtomicShell* shell = GetAtomicShell(Z, as);
	G4double bindingEnergy = shell->BindingEnergy();

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
	      stot -= mfp*std::log(G4UniformRand());
	      if( stot > 1.0 || eLossMax < bindingEnergy) { break; }
	      // sample deexcitation
	      vdyn.clear();
	      GenerateParticles(&vdyn, shell, Z, gCut, eCut); 
	      G4int nsec = vdyn.size();
	      if(nsec > 0) {
		G4ThreeVector r = prePos  + stot*delta;
		G4double time   = preTime + stot*dt;
		for(G4int j=0; j<nsec; ++j) {
		  G4DynamicParticle* dp = vdyn[j];
		  G4double e = dp->GetKineticEnergy();

		  // save new secondary if there is enough energy
		  if(eLossMax >= e) {
		    eLossMax -= e;		    
		    G4Track* t = new G4Track(dp, time, r);
		    tracks.push_back(t);
		  } else {
		    delete dp;
		  }
		}
	      }
	    } while (stot < 1.0);
	  }
	}
      }
    } 
  }
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
