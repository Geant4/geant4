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
// $Id: G4VIStore.hh,v 1.7 2002-09-02 13:25:26 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VIStore
//
// Class description:
//
// An interface of an "importance store" used by importance sampling.
// It defines how a importance value together with a "cell" 
// (a G4VPhysicalVolume and a replica number) has to be added
// to the "importance store" and how a importance value can be derived 
// from the "importance store". 
// 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VIStore_hh
#define G4VIStore_hh G4VIStore_hh

#include "globals.hh"

class G4GeometryCell;
class G4VPhysicalVolume;

class  G4VIStore
{

public:  // with description

  virtual  ~G4VIStore() {}

  virtual G4double GetImportance(const G4GeometryCell &gCell) const = 0;
    // derive a importance value of a "cell" addresed by a G4GeometryCell
    // from the store.

  virtual G4bool IsKnown(const G4GeometryCell &gCell) const = 0;
    // returns true if the gCell is in the store, else false 

  virtual const G4VPhysicalVolume &GetWorldVolume() const = 0;
};

#endif
