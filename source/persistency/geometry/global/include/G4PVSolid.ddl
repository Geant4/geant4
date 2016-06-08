// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PVSolid.ddl,v 1.3 1999/12/15 14:51:22 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// 
// class G4PVSolid
//
// History:
// 19.06.98 A.Kimura Converted from G4VSolid.hh
// 20.06.98 A.Kimura Clean up

#ifndef G4PVSOLID_DDL
#define G4PVSOLID_DDL

#include "globals.hh"
#include "geomdefs.hh"
#include "G4PersistentTypes.hh"
#include "G4PersistentSchema.hh"

#include "G4VSolid.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class G4PVSolid;

class G4PVSolid : public HepPersObj {
public:
    G4PVSolid();
    G4PVSolid(const G4String& name);
    virtual ~G4PVSolid();

    virtual G4VSolid* MakeTransientObject() const = 0;

    G4bool operator==( const G4PVSolid& s) const
    {
	return (this==&s) ? true : false;
    }

    G4String GetName() const;
    void SetName(const G4String& name);
                                   
    virtual G4GeometryType  GetEntityType() const = 0;

private:

    G4PString fshapeName;	 // Name
};

#endif
