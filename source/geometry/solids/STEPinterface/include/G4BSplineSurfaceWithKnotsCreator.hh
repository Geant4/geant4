// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BSplineSurfaceWithKnotsCreator.hh,v 1.3 2000-11-09 16:35:43 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BSplineSurfaceWithKnotsCreator
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
#ifndef G4BSPLINESURFACEWITHKNOTSCREATOR_HH
#define G4BSPLINESURFACEWITHKNOTSCREATOR_HH

#include "G4BSplineSurfaceCreator.hh"

class G4BSplineSurfaceWithKnotsCreator: public G4BSplineSurfaceCreator
{
  public:

  // Constructor & destructor

    G4BSplineSurfaceWithKnotsCreator();
    ~G4BSplineSurfaceWithKnotsCreator();

  // Member functions

    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void*);
    const char* Name() const { return "B_Spline_Surface_With_Knots"; }
    static G4BSplineSurfaceWithKnotsCreator GetInstance() { return csc; }

  // Members

  private:

    static G4BSplineSurfaceWithKnotsCreator csc;
};

#endif
