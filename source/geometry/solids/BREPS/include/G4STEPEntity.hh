// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4STEPEntity.hh,v 1.3 2000-08-28 15:00:34 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4STEPEntity
//
// Class description:
// 
// Utility abstract class for STEP entity.

// Authors: J.Sulkimo, P.Urban.
// ----------------------------------------------------------------------
#ifndef __G4STEPENTITY
#define __G4STEPENTITY

#include "globals.hh"
#include "G4OrderedTable.hh"

class G4STEPEntity
{

public:

  virtual G4String GetEntityType() const = 0;

};

#endif
