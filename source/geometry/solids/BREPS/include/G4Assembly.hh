// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Assembly.hh,v 1.6 2001-02-05 18:45:32 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include "G4PlacedSolid.hh"
#include "G4BREPSolid.hh"

typedef G4RWTPtrOrderedVector<G4PlacedSolid> G4PlacedVector;  

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
