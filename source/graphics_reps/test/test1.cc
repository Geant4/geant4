// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test1.cc,v 1.1 1999-01-08 16:32:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  19th July 1996

// Simply tests linkability of kernel without visualization.

#include "G4ios.hh"
#include <math.h>

#include "BuildCalorimeter.hh"
#include "G4VPhysicalVolume.hh"
#include "G4GeometryManager.hh"

#include "g4templates.hh"

main () {

  // Build some geometry and get relevant pointers.
  G4VPhysicalVolume *pWorld = BuildCalorimeter ();
  G4GeometryManager::GetInstance () -> CloseGeometry ();

  return 0;
}
