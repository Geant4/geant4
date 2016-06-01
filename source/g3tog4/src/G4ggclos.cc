// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ggclos.cc,v 2.1 1998/07/13 16:50:30 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
#include "G4GeometryManager.hh"
#include "G3G4Interface.hh"

void PG4ggclos()
{
  G4ggclos();
}

void G4ggclos()
{
        // close out the geometry
    
    G4cout << "Closing geometry..." << endl;
  G4bool optimise=false;
  G4GeometryManager::GetInstance()->CloseGeometry(optimise);
  G4cout << "Geometry closed." << endl;
}
