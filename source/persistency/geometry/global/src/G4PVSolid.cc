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
// $Id: G4PVSolid.cc,v 1.3 2001/07/11 10:02:18 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// class G4PVSolid
//
// Implementation for persistent solid base class
//
//
// History:
//  19.06.98 A.Kimura Converted G4VSolid.cc

#include "G4PVSolid.hh"

// Constructor (default)
G4PVSolid::G4PVSolid()
 : fshapeName("____________Unknown_____________") 
{
}

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
