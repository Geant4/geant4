// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ggclos.cc,v 1.1 1999-01-07 16:06:48 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
