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
// $Id: G4VWeightWindowStore.hh,v 1.2 2003/08/19 15:44:57 dressel Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// ----------------------------------------------------------------------
// Class G4VWeightWindowStore
//
// Class description:
//
// Interface class for a weight window store. It defines how the lower 
// weight window bound can be obtained from a weight window store.
// 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VWeightWindowStore_hh
#define G4VWeightWindowStore_hh G4VWeightWindowStore_hh

#include "globals.hh"

class G4GeometryCell;
class G4VPhysicalVolume;

class  G4VWeightWindowStore
{

public:  // with description

  G4VWeightWindowStore();
  virtual  ~G4VWeightWindowStore();

  virtual G4double GetLowerWeitgh(const G4GeometryCell &gCell, 
			 G4double partEnergy) const = 0;
    // derive a lower weight bound value of a "cell" addresed by a 
    // G4GeometryCell and the coresponding energy from the store.

  virtual G4bool IsKnown(const G4GeometryCell &gCell) const = 0;
    // returns true if the gCell is in the store, else false 


  virtual const G4VPhysicalVolume &GetWorldVolume() const = 0;
    // return a reference to the wolrd volume of the 
    // geometry
};

#endif
