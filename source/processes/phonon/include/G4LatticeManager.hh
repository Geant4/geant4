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
// $Id$
//
#ifndef G4LatticeManager_h
#define G4LatticeManager_h 1


#include "G4ThreeVector.hh"
#include <map>

class G4LatticeLogical;
class G4LatticePhysical;
class G4Material;
class G4VPhysicalVolume;


class G4LatticeManager {
private:
  static G4ThreadLocal G4LatticeManager* fLM;		// Singleton

public:
  static G4LatticeManager* GetLatticeManager(); 

  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }

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
  G4int verboseLevel;		// Allow users to enable diagnostic messages

  typedef std::map<G4VPhysicalVolume*, G4LatticePhysical*> LatticeVolMap;
  typedef std::map<G4Material*, G4LatticeLogical*> LatticeMatMap;

  LatticeMatMap fLLatticeList;
  int fTotalLLattices;		// == fLLatticeList.size(), for convenience

  LatticeVolMap fPLatticeList; 
  int fTotalPLattices;		// == fPLatticeList.size(), for convenience

private:
  G4LatticeManager();
  virtual ~G4LatticeManager();
};

#endif
