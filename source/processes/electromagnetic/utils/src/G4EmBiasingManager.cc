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
// $Id: G4EmBiasingManager.cc,v 1.88 2010-08-17 17:36:59 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4EmBiasingManager
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 28.07.2011
//
// Modifications:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmBiasingManager.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmBiasingManager::G4EmBiasingManager() 
  : nRegions(0),currentStepLimit(0.0),startTracking(true)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmBiasingManager::~G4EmBiasingManager()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmBiasingManager::Initialise()
{
  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  idxCouple.resize(numOfCouples, -1);

  // Deexcitation
  if (nRegions>0) {

    for (size_t j=0; j<numOfCouples; ++j) {
      const G4MaterialCutsCouple* couple =
        theCoupleTable->GetMaterialCutsCouple(j);
      const G4ProductionCuts* pcuts = couple->GetProductionCuts();
      for(G4int i=0; i<nRegions; ++i) {
	if(forcedRegions[i]) {
	  if(pcuts == forcedRegions[i]->GetProductionCuts()) { 
	    idxCouple[j] = i;
	    break; 
	  }
	}
      }
    }
  }
  if (nRegions > 0) {
    G4cout << " Forced Interaction is activated for G4Regions: " << G4endl;
    for (G4int i=0; i<nRegions; ++i) {
      const G4Region* r = forcedRegions[i];
      if(r) { G4cout << "           " << r->GetName() << G4endl; }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmBiasingManager::ActivateForcedInteraction(G4double val, 
						   const G4String& rname)
{
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  G4String name = rname;
  if(name == "" || name == "world" || name == "World") {
    name = "DefaultRegionForTheWorld";
  }
  const G4Region* reg = regionStore->GetRegion(name, false);
  if(!reg) { 
    G4cout << "### G4EmBiasingManager::ActivateDeexcitation WARNING: "
	   << " G4Region <"
	   << rname << "> is unknown" << G4endl;
    return; 
  }

  // the region is in the list
  if (nRegions) {
    for (G4int i=0; i<nRegions; ++i) {
      if (reg == forcedRegions[i]) {
	lengthForRegion[i] = val; 
        return;
      }
    }
  }
  if(val < 0.0) { 
    G4cout << "### G4EmBiasingManager::ActivateDeexcitation WARNING: "
	   << val << " < 0.0, so no activation for the G4Region <"
	   << rname << ">" << G4endl;
    return; 
  }

  // new region 
  forcedRegions.push_back(reg);
  lengthForRegion.push_back(val);
  ++nRegions;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmBiasingManager::GetStepLimit(G4int coupleIdx, 
					  G4double previousStep)
{
  if(startTracking) {
    startTracking = false;
    G4int i = idxCouple[coupleIdx];
    if(i < 0) {
      currentStepLimit = DBL_MAX;
    } else {
      currentStepLimit = lengthForRegion[i];
      if(currentStepLimit > 0.0) { currentStepLimit *= G4UniformRand(); }
    }
  } else {
    currentStepLimit -= previousStep;
  }
  return currentStepLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
