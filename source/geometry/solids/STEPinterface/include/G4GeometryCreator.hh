// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GeometryCreator.hh,v 1.5 2000-11-10 17:44:45 gcosmo Exp $
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
  
    G4GeometryCreator();
    virtual ~G4GeometryCreator();

  // Member functions
  
    virtual G4bool operator==(const G4GeometryCreator&);

    virtual void CreateG4Geometry(STEPentity&)=0;
    virtual void CreateSTEPGeometry(void* =0)=0;

    virtual void* GetCreatedObject();
  
    virtual const char* Name() const=0;
    virtual STEPattribute* GetNamedAttribute(const G4String&, STEPentity&);
    virtual STEPentity* GetNamedEntity(const G4String&, STEPentity&);  
    InstMgr* GetInstanceManager() const;
    G4int GetNextId() { objectId+=10; return objectId; }

  protected:
    
    static G4int objectId;
    static InstMgr instanceManager;
    void* createdObject;

};

#endif
