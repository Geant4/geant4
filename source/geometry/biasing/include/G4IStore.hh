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
// $Id: G4IStore.hh,v 1.7 2002-09-02 13:25:25 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4IStore
//
// Class description:
//
// An implementation of a "importance store" with the interface
// G4VIStore. See description in G4VIStore.
// This implementation uses G4GeometryCellImportance as the container
// to store the "cells" together with the importance values.
// Giving a cell the importance 0 is allowed as
// a flagging that no biasing should happen between this
// cell and it's neighbors
// If a cell is not known by the importance store no biasing
// should be applied between this cell and it's nighbors.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4IStore_hh
#define G4IStore_hh G4IStore_hh 

#include "G4VIStore.hh"
#include "G4GeometryCellImportance.hh"

class G4IStore : public G4VIStore
{

public:  // with description

  G4IStore(const G4VPhysicalVolume &worldvolume);
    // initialise the importance store for the given geometry

  ~G4IStore();
    // destruct

  void AddImportanceGeometryCell(G4double importance,
			   const G4GeometryCell &gCell);
  void AddImportanceGeometryCell(G4double importance,
			   const G4VPhysicalVolume &,
			   G4int aRepNum = 0);
    // Add a "cell" together with a importance value to the store.

  void ChangeImportance(G4double importance,
			const G4GeometryCell &gCell);
  void ChangeImportance(G4double importance,
			const G4VPhysicalVolume &,
			G4int aRepNum = 0);
    // Change a importance value of a "cell".

  G4double GetImportance(const G4VPhysicalVolume &,
			 G4int aRepNum = 0) const ;
  G4double GetImportance(const G4GeometryCell &gCell) const;
    // derive a importance value of a "cell" addresed by a G4GeometryCell
    // from the store.

  G4bool IsKnown(const G4GeometryCell &gCell) const;
    // returns true if the gCell is in the store, else false 


  const G4VPhysicalVolume &GetWorldVolume() const;
    // return a reference to the wolrd volume of the 
    // "importance" geometry
  
private:

  G4bool IsInWorld(const G4VPhysicalVolume &) const;
  void SetInternalIterator(const G4GeometryCell &gCell) const;
  void Error(const G4String &m) const;
  
private:
 
  const G4VPhysicalVolume &fWorldVolume;
  G4GeometryCellImportance fGeometryCelli;

  mutable G4GeometryCellImportance::const_iterator fCurrentIterator;
};

#endif
