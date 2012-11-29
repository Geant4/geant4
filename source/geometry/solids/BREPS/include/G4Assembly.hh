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
// $Id$
//
// ----------------------------------------------------------------------
// Class G4Assembly
//
// Class description:
//   
// Object providing access to the final assembled geometry as read, for
// instance, from a STEP file description. It owns a vector of placed
// solids which must be initialised through the SetPlacedVector()
// method. Currently, it simply provides a way to retrieve the pointer
// to each placed solid, and to determine the number of solids currently
// placed.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------

#ifndef G4ASSEMBLY_HH
#define G4ASSEMBLY_HH

#include <vector>
#include "G4PlacedSolid.hh"
#include "G4BREPSolid.hh"

typedef std::vector<G4PlacedSolid*> G4PlacedVector;  

class G4Assembly
{

public: // with description

  G4Assembly();
  ~G4Assembly();
    // Constructor & destructor

  void SetPlacedVector(G4PlacedVector&);
  inline G4PlacedSolid* GetPlacedSolid(G4int solidNumber) const;
  inline G4int GetNumberOfSolids() const;

private:

  G4Assembly(const G4Assembly&);
  G4Assembly& operator=(const G4Assembly&);
    // Private copy constructor and assignment operator.

private:  

  G4int           numberOfSolids;
  G4PlacedVector  placedVec;

};

#include "G4Assembly.icc"

#endif
