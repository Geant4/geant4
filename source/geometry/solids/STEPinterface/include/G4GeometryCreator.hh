// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GeometryCreator.hh,v 1.3 2000-01-21 13:45:26 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4GeometryCreator
//
// Class description:
//
//

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------
#ifndef G4GEOMETRYCREATOR_HH
#define G4GEOMETRYCREATOR_HH

#include <schema.h>
#include "globals.hh"
#include "G4PlacedSolid.hh"
#include "G4Surface.hh"
#include "G4BREPSolid.hh"

typedef G4RWTPtrOrderedVector<G4PlacedSolid> G4PlacedSolidVector;
typedef G4RWTPtrOrderedVector<G4Surface> G4SurfaceVector;
typedef G4RWTPtrOrderedVector<G4BREPSolid> G4SolidVector;

class G4GeometryCreator
{
  
  public:

  // Constructor & destructor
  
    G4GeometryCreator() {;}
    virtual ~G4GeometryCreator() {;}

  // Member functions
  
    virtual void CreateG4Geometry(STEPentity&)=0;
    virtual void CreateSTEPGeometry(void* =0)=0;

    virtual void* GetCreatedObject() { return createdObject; }
  
    virtual G4String Name()=0;
    virtual G4bool operator==(const G4GeometryCreator&) { return 0; }
    virtual STEPattribute* GetNamedAttribute(G4String&,STEPentity&);
    virtual STEPentity* GetNamedEntity(G4String&,STEPentity&);  
    G4int GetNextId() { objectId+=10; return objectId; }

  // Members
    
    static G4int objectId;
    static InstMgr instanceManager;
    void* createdObject;

};

#endif
