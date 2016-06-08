// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PVSolid.cc,v 1.1 1999/01/07 16:10:49 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// class G4PVSolid
//
// Implementation for persistent solid base class
//
//
// History:
//  19.06.98 A.Kimura Converted G4VSolid.cc

#include "G4PVSolid.hh"

// Constructor
//  - Copies name
G4PVSolid::G4PVSolid(const G4String& name)
 : fshapeName(name) 
{
}

// Destructor (virtual)
// - Remove ourselves from solid Store
G4PVSolid::~G4PVSolid()
{
}
	
// Returns name by value
G4String G4PVSolid::GetName() const
{
    return (const char*) fshapeName;
}

void G4PVSolid::SetName(const G4String& name)
{
    fshapeName = name;
}	
