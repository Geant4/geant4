// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RationalBSplineSurfaceCreator.hh,v 1.2 2000-01-21 13:45:29 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4RationalBSplineSurfaceCreator
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
#ifndef G4RATIONALBSPLINESURFACECREATOR_HH
#define G4RATIONALBSPLINESURFACECREATOR_HH

#include "G4GeometryCreator.hh"
#include "G4BSplineSurfaceCreator.hh"

class G4RationalBSplineSurfaceCreator: public G4BSplineSurfaceCreator 
{
  public:

  // Constructor & destructor

    G4RationalBSplineSurfaceCreator();
    ~G4RationalBSplineSurfaceCreator();

  // Member functions

    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void*);
    G4String Name() { return "Rational_B_Spline_Surface"; }
    static G4RationalBSplineSurfaceCreator GetInstance() { return csc; }
    
  // Members

  private:

    static G4RationalBSplineSurfaceCreator csc;
};

#endif
