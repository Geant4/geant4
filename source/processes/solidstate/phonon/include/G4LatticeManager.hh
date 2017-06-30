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
/// \file materials/include/G4LatticeManager.hh
/// \brief Definition of the G4LatticeManager class
//
// $Id: G4LatticeManager.hh 84149 2014-10-08 18:04:16Z mkelsey $
//
// 20131113  Add registry to carry unique lattice pointers, for EOJ deletion
// 20131115  Drop lattice counters, not used anywhere
// 20141008  Change to global singleton; must be shared across worker threads

#ifndef G4LatticeManager_h
#define G4LatticeManager_h 1

#include "G4ThreeVector.hh"
#include <map>
#include <set>

class G4LatticeLogical;
class G4LatticePhysical;
class G4Material;
class G4VPhysicalVolume;


class G4LatticeManager {
private:
  static G4LatticeManager* fLM;		// Global, shared singleton

public:
  static G4LatticeManager* GetLatticeManager(); 

  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }

  void Reset();		// Remove and delete all registered lattices

  // Users may register physical or logical lattices with volumes
  G4bool RegisterLattice(G4VPhysicalVolume*, G4LatticePhysical*);
  G4bool RegisterLattice(G4VPhysicalVolume*, G4LatticeLogical*);

  // Logical lattices are associated with materials
  G4bool RegisterLattice(G4Material*, G4LatticeLogical*);

  // Logical lattices may be read from <latDir>/config.txt data file
  G4LatticeLogical* LoadLattice(G4Material*, const G4String& latDir);
  G4LatticeLogical* GetLattice(G4Material*) const;
  G4bool HasLattice(G4Material*) const;

  // Combine loading and registration (Material extracted from volume)
  G4LatticePhysical* LoadLattice(G4VPhysicalVolume*, const G4String& latDir);

  // NOTE:  Passing Vol==0 will return the default lattice
  G4LatticePhysical* GetLattice(G4VPhysicalVolume*) const;
  G4bool HasLattice(G4VPhysicalVolume*) const;

  G4double MapKtoV(G4VPhysicalVolume*, G4int, const G4ThreeVector &) const;

  G4ThreeVector MapKtoVDir(G4VPhysicalVolume*, G4int,
			   const G4ThreeVector&) const;

protected:
  void Clear();		// Remove entries from lookup tables w/o deletion

protected:
  G4int verboseLevel;		// Allow users to enable diagnostic messages

  typedef std::map<G4Material*, G4LatticeLogical*> LatticeMatMap;
  typedef std::set<G4LatticeLogical*> LatticeLogReg;

  LatticeLogReg fLLattices;	// Registry of unique lattice pointers
  LatticeMatMap fLLatticeList;

  typedef std::map<G4VPhysicalVolume*, G4LatticePhysical*> LatticeVolMap;
  typedef std::set<G4LatticePhysical*> LatticePhyReg;

  LatticePhyReg fPLattices;	// Registry of unique lattice pointers
  LatticeVolMap fPLatticeList; 

private:
  G4LatticeManager();
  virtual ~G4LatticeManager();
};

#endif
