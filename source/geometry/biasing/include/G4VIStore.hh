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
//
// $Id: G4VIStore.hh 66356 2012-12-18 09:02:32Z gcosmo $
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

  G4VIStore();
  virtual  ~G4VIStore();

  virtual G4double GetImportance(const G4GeometryCell &gCell) const = 0;
    // derive a importance value of a "cell" addresed by a G4GeometryCell
    // from the store.

  virtual G4bool IsKnown(const G4GeometryCell &gCell) const = 0;
    // returns true if the gCell is in the store, else false 

  virtual const G4VPhysicalVolume &GetWorldVolume() const = 0;
};

#endif
