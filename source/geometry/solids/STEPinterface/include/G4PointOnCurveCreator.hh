// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PointOnCurveCreator.hh,v 1.2 2000-01-21 13:45:29 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4PointOnCurveCreator
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
#ifndef G4POINTONCURVECREATOR_HH
#define G4POINTONCURVECREATOR_HH

#include "G4GeometryCreator.hh"

class G4PointOnCurveCreator: private G4GeometryCreator 
{
  public:

  // Constructor & destructor

    G4PointOnCurveCreator();
    ~G4PointOnCurveCreator();

  // Member functions

    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void* G4obj);
    G4String Name() { return "Point_On_Curve"; }
    static G4PointOnCurveCreator GetInstance() { return csc; }

  // Members

  private:

    static G4PointOnCurveCreator csc;
};

#endif
