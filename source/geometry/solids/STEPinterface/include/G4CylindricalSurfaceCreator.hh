// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4CylindricalSurfaceCreator.hh,v 1.3 2000-11-09 16:35:44 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4CylindricalSurfaceCreator
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
#ifndef G4CYLINDRICALSURFACECREATOR_HH
#define G4CYLINDRICALSURFACECREATOR_HH

#include "G4GeometryCreator.hh"

class G4CylindricalSurfaceCreator: private G4GeometryCreator 
{
  public:

  // Constructor & destructor

    G4CylindricalSurfaceCreator();
    ~G4CylindricalSurfaceCreator();

  // Member functions

    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void* G4obj);
    const char* Name() const { return "Cylindrical_Surface"; }
    static G4CylindricalSurfaceCreator GetInstance() { return csc; }

  // Members
  
  private:

    static G4CylindricalSurfaceCreator csc;
};

#endif
