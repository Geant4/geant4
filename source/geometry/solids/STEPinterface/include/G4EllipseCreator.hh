// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EllipseCreator.hh,v 1.3 2000-11-09 16:35:45 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4EllipseCreator
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
#ifndef G4ELLIPSECREATOR_HH
#define G4ELLIPSECREATOR_HH

#include "G4GeometryCreator.hh"

class G4EllipseCreator: private G4GeometryCreator 
{
  public:
  
  // Constructor & destructor

    G4EllipseCreator();
    ~G4EllipseCreator();

  // Member functions

    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void* G4obj);
    const char* Name() const { return "Ellipse"; }
    static G4EllipseCreator GetInstance() { return csc; }

  //Members

  private:

    static G4EllipseCreator csc;
};

#endif
