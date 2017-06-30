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
/// \file materials/src/G4LatticeManager.cc
/// \brief Implementation of the G4LatticeManager class
//
// $Id: G4LatticeManager.cc 84149 2014-10-08 18:04:16Z mkelsey $
//
// 20131113  Delete lattices in (new) registry, not in lookup maps
// 20141008  Change to global singleton; must be shared across worker threads

#include "G4LatticeManager.hh"
#include "G4AutoLock.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticePhysical.hh"
#include "G4LatticeReader.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include <fstream>

G4LatticeManager* G4LatticeManager::fLM = 0;	// Global (shared) singleton

namespace {
  G4Mutex latMutex = G4MUTEX_INITIALIZER;	// For thread protection
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LatticeManager::G4LatticeManager() : verboseLevel(0) {
  Clear();
}

G4LatticeManager::~G4LatticeManager() {
  Reset();			// Deletes all lattices
}

// Delete all registered lattices and clear entries from lookup tables

void G4LatticeManager::Reset() {
  for (LatticeLogReg::iterator lm=fLLattices.begin();
       lm != fLLattices.end(); ++lm) {
    delete (*lm);
  }

  for (LatticePhyReg::iterator pm=fPLattices.begin();
       pm != fPLattices.end(); ++pm) {
    delete (*pm);
  }

  Clear();
}

// Remove entries without deletion (for begin-job and end-job initializing)

void G4LatticeManager::Clear() {
  fPLatticeList.clear();
  fPLattices.clear();

  fLLatticeList.clear();
  fLLattices.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LatticeManager* G4LatticeManager::GetLatticeManager() {
  // if no lattice manager exists, create one.
  G4AutoLock latLock(&latMutex);	// Protect before changing pointer
  if (!fLM) fLM = new G4LatticeManager();
  latLock.unlock();

  return fLM;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Associate logical lattice with material

G4bool G4LatticeManager::RegisterLattice(G4Material* Mat,
					 G4LatticeLogical* Lat) {
  if (!Mat || !Lat) return false;	// Don't register null pointers

  G4AutoLock latLock(&latMutex);	// Protect before changing registry
  fLLattices.insert(Lat);		// Take ownership in registry
  fLLatticeList[Mat] = Lat;
  latLock.unlock();

  if (verboseLevel) {
    G4cout << "G4LatticeManager::RegisterLattice: "
	   << " Total number of logical lattices: " << fLLatticeList.size()
	   << " (" << fLLattices.size() << " unique)" << G4endl;
  }

  return true; 
}

// Construct logical lattice for material from config file

G4LatticeLogical* G4LatticeManager::LoadLattice(G4Material* Mat,
						const G4String& latDir) {
  if (verboseLevel) {
    G4cout << "G4LatticeManager::LoadLattice material " << Mat->GetName()
	   << " " << latDir << G4endl;
  }
		      
  G4LatticeReader latReader(verboseLevel);
  G4LatticeLogical* newLat = latReader.MakeLattice(latDir+"/config.txt");
  if (verboseLevel>1) G4cout << " Created newLat " << newLat << G4endl;

  if (newLat) RegisterLattice(Mat, newLat);
  else {
    G4cerr << "ERROR creating " << latDir << " lattice for material "
	   << Mat->GetName() << G4endl;
  }

  return newLat;
}

// Combine loading and registration (Material extracted from volume)

G4LatticePhysical* G4LatticeManager::LoadLattice(G4VPhysicalVolume* Vol,
						const G4String& latDir) {
  if (verboseLevel) {
    G4cout << "G4LatticeManager::LoadLattice volume " << Vol->GetName()
	   << " " << latDir << G4endl;
  }
		      
  G4Material* theMat = Vol->GetLogicalVolume()->GetMaterial();

  // Create and register the logical lattice, then the physical lattice
  G4LatticeLogical* lLattice = LoadLattice(theMat, latDir);
  if (!lLattice) return 0;

  G4LatticePhysical* pLattice =
    new G4LatticePhysical(lLattice, Vol->GetFrameRotation());
  if (pLattice) RegisterLattice(Vol, pLattice);

  if (verboseLevel>1) G4cout << " Created pLattice " << pLattice << G4endl;

  return pLattice;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Associate physical (oriented) lattice with physical volume

G4bool G4LatticeManager::RegisterLattice(G4VPhysicalVolume* Vol,
					 G4LatticePhysical* Lat) {
  if (!Vol || !Lat) return false;	// Don't register null pointers

  G4AutoLock latLock(&latMutex);	// Protect before changing registry

  // SPECIAL: Register first lattice with a null volume to act as default
  if (fPLatticeList.empty()) fPLatticeList[0] = Lat;

  fPLattices.insert(Lat);
  fPLatticeList[Vol] = Lat;

  latLock.unlock();

  if (verboseLevel) {
    G4cout << "G4LatticeManager::RegisterLattice: "
	   << " Total number of physical lattices: " << fPLatticeList.size()-1
	   << " (" << fPLattices.size() << " unique)" << G4endl;
  }

  return true; 
}

G4bool G4LatticeManager::RegisterLattice(G4VPhysicalVolume* Vol,
					 G4LatticeLogical* LLat) {
  if (!Vol || !LLat) return false;	// Don't register null pointers

  // Make sure logical lattice is registered for material
  RegisterLattice(Vol->GetLogicalVolume()->GetMaterial(), LLat);

  // Create and register new physical lattice to go with volume
  return RegisterLattice(Vol, new G4LatticePhysical(LLat, Vol->GetFrameRotation()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Returns a pointer to the LatticeLogical associated with material

G4LatticeLogical* G4LatticeManager::GetLattice(G4Material* Mat) const {
  LatticeMatMap::const_iterator latFind = fLLatticeList.find(Mat);
  if (latFind != fLLatticeList.end()) {
    if (verboseLevel)
      G4cout << "G4LatticeManager::GetLattice found " << latFind->second
	     << " for " << (Mat?Mat->GetName():"NULL") << "." << G4endl;
    return latFind->second;
  }

  if (verboseLevel) 
    G4cerr << "G4LatticeManager:: Found no matching lattices for "
	   << (Mat?Mat->GetName():"NULL") << "." << G4endl;

  return 0;			// No lattice associated with volume
}

// Returns a pointer to the LatticePhysical associated with volume
// NOTE:  Passing Vol==0 will return the default lattice

G4LatticePhysical* G4LatticeManager::GetLattice(G4VPhysicalVolume* Vol) const {
  LatticeVolMap::const_iterator latFind = fPLatticeList.find(Vol);
  if (latFind != fPLatticeList.end()) {
    if (verboseLevel)
      G4cout << "G4LatticeManager::GetLattice found " << latFind->second
	     << " for " << (Vol?Vol->GetName():"default") << "." << G4endl;
    return latFind->second;
  }

  if (verboseLevel) 
    G4cerr << "G4LatticeManager::GetLattice found no matching lattices for "
	   << (Vol?Vol->GetName():"default") << "." << G4endl;

  return 0;			// No lattice associated with volume
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Return true if volume Vol has a physical lattice

G4bool G4LatticeManager::HasLattice(G4VPhysicalVolume* Vol) const {
  return (fPLatticeList.find(Vol) != fPLatticeList.end());
}

// Return true if material Mat has a logical lattice

G4bool G4LatticeManager::HasLattice(G4Material* Mat) const {
  return (fLLatticeList.find(Mat) != fLLatticeList.end());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//Given the phonon wave vector k, phonon physical volume Vol 
//and polarizationState(0=LON, 1=FT, 2=ST), 
//returns phonon velocity in m/s

G4double G4LatticeManager::MapKtoV(G4VPhysicalVolume* Vol,
				 G4int polarizationState,
				 const G4ThreeVector & k) const {
  G4LatticePhysical* theLattice = GetLattice(Vol);
  if (verboseLevel)
    G4cout << "G4LatticeManager::MapKtoV using lattice " << theLattice
	   << G4endl;

  // If no lattice available, use generic "speed of sound"
  return theLattice ? theLattice->MapKtoV(polarizationState, k) : 300.*m/s;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Given the phonon wave vector k, phonon physical volume Vol 
// and polarizationState(0=LON, 1=FT, 2=ST), 
// returns phonon propagation direction as dimensionless unit vector

G4ThreeVector G4LatticeManager::MapKtoVDir(G4VPhysicalVolume* Vol,
					   G4int polarizationState,
					   const G4ThreeVector & k) const {
  G4LatticePhysical* theLattice = GetLattice(Vol);
  if (verboseLevel)
    G4cout << "G4LatticeManager::MapKtoVDir using lattice " << theLattice
	   << G4endl;

  // If no lattice available, propagate along input wavevector
  return theLattice ? theLattice->MapKtoVDir(polarizationState, k) : k.unit();
}
