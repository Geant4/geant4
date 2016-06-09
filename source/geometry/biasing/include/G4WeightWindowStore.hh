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
// $Id: G4WeightWindowStore.hh,v 1.2 2003/08/19 15:44:57 dressel Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// ----------------------------------------------------------------------
// Class G4WeightWindowStore
//
// Class description:
//
// Implementation of a weight window store according to the
// G4VWeightWindowStore interface.
// See G4VWeightWindowStore.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4WeightWindowStore_hh
#define G4WeightWindowStore_hh G4WeightWindowStore_hh 

#include "G4VWeightWindowStore.hh"
#include "G4GeometryCellWeight.hh"
#include <set>
#include <vector>

class G4WeightWindowStore: public G4VWeightWindowStore
{

public:  // with description

  explicit G4WeightWindowStore(const G4VPhysicalVolume &worldvolume);
    // initialise the weight window store for the given geometry

  virtual ~G4WeightWindowStore();
    // destruct

  virtual G4double GetLowerWeitgh(const G4GeometryCell &gCell, 
			 G4double partEnergy) const;
    // derive a lower weight bound value of a "cell" addresed by a 
    // G4GeometryCell and the coresponding energy from the store.

  virtual G4bool IsKnown(const G4GeometryCell &gCell) const;
    // returns true if the gCell is in the store, else false 


  virtual const G4VPhysicalVolume &GetWorldVolume() const;
    // return a reference to the wolrd volume of the 
    // geometry

  void AddLowerWeights(const G4GeometryCell &gCell,
		       const std::vector<G4double> &lowerWeights);
    // add lower weights. Only if general upper energy bounds have 
    // been set
 
  void AddUpperEboundLowerWeightPairs(const G4GeometryCell &gCell,
				      const G4UpperEnergyToLowerWeightMap&
				      enWeMap);
    // set upper energy - lower weight pairs for a cell

  void SetGeneralUpperEnergyBounds(const std::set<G4double, std::less<G4double> > &enBounds);
  
private:

  G4bool IsInWorld(const G4VPhysicalVolume &) const;
  void Error(const G4String &m) const;
  void SetInternalIterator(const G4GeometryCell &gCell) const;
  
private:

  const G4VPhysicalVolume &fWorldVolume;  
  std::set<G4double, std::less<G4double> > fGeneralUpperEnergyBounds;
  G4GeometryCellWeight fCellToUpEnBoundLoWePairsMap;
  mutable G4GeometryCellWeight::const_iterator fCurrentIterator;
  
};

#endif
