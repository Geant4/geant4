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
// $Id: G4GeometryCreator.hh,v 1.7 2001-07-11 10:00:05 gunter Exp $
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
#include "g4std/vector"
#include "globals.hh"
#include "G4PlacedSolid.hh"
#include "G4Surface.hh"
#include "G4BREPSolid.hh"

typedef G4std::vector<G4PlacedSolid*> G4PlacedSolidVector;
typedef G4std::vector<G4Surface*> G4SurfaceVector;
typedef G4std::vector<G4BREPSolid*> G4SolidVector;

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
