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
// $Id: G4VIStore.hh,v 1.3 2002-04-10 13:13:07 dressel Exp $
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

class G4PTouchableKey;
class G4VPhysicalVolume;

class  G4VIStore
{

public:  // with description

  G4VIStore(G4VPhysicalVolume &worldVolume){}
  virtual  ~G4VIStore(){}

  virtual void AddImportanceRegion(G4double importance,
				   const G4VPhysicalVolume &,
				   G4int aRepNum = 0) = 0;
    // Add a "cell" together with a importance value to the store.

  virtual void ChangeImportance(G4double importance,
				const G4VPhysicalVolume &,
				G4int aRepNum = 0) = 0;
    // Change a importance value of a "cell".

  virtual G4double GetImportance(const G4VPhysicalVolume &,
				 G4int aRepNum = 0) const = 0;
    // derive the importance value of a "cell" from the store.

  virtual G4double GetImportance(const G4PTouchableKey &ptk) const = 0;
    // derive a importance value of a "cell" addresed by a G4PTouchableKey
    // from the store.

  virtual G4VPhysicalVolume &GetWorldVolume() = 0;
};

#endif
