//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4IStore.hh,v 1.3 2002-04-10 13:13:06 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4IStore
//
// Class description:
//
// An implementation of a "importance store" with the interface
// G4VIStore. See description in G4VIStore.
// This implementation uses G4PtkImportance as the container
// to store the "cells" together with the importance values.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4IStore_hh
#define G4IStore_hh G4IStore_hh 

#include "G4VIStore.hh"
#include "G4PtkImportance.hh"

class G4IStore : public G4VIStore
{

public:  // with description

  G4IStore(G4VPhysicalVolume &worldvolume);
    // initialise the importance store for the given geometry

  ~G4IStore();
    // destruct

  void AddImportanceRegion(G4double importance,
			   const G4VPhysicalVolume &,
			   G4int aRepNum = 0);
    // Add a "cell" together with a importance value to the store.

  void ChangeImportance(G4double importance,
			const G4VPhysicalVolume &,
			G4int aRepNum = 0);
    // Change a importance value of a "cell".

  G4double GetImportance(const G4VPhysicalVolume &,
			 G4int aRepNum = 0) const ;
    // derive the importance value of a "cell" from the store.

  G4double GetImportance(const G4PTouchableKey &ptk) const;
    // derive a importance value of a "cell" addresed by a G4PTouchableKey
    // from the store.

  G4VPhysicalVolume &GetWorldVolume();
    // return a reference to the wolrd volume of the 
    // "importance" geometry
  
private:

  G4bool IsInWorld(const G4VPhysicalVolume &) const;
  void SetInternalIterator(const G4VPhysicalVolume &,
			   G4int aRepNum) const;
  void Error(const G4String &m) const;

private:
 
  G4VPhysicalVolume &fWorldVolume;
  G4PtkImportance fPtki;

  mutable G4PtkImportance::const_iterator fCurrentIterator;
};

#endif
