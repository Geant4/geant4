// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AdvancedBrepShapeRepresentationCreator.hh,v 1.3 2000-11-09 16:35:41 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4AdvancedBrepShapeRepresentationCreator
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
#ifndef G4ADVANCEDBREPSHAPEREPRESENTATIONCREATOR_HH
#define G4ADVANCEDBREPSHAPEREPRESENTATIONCREATOR_HH

#include "G4GeometryCreator.hh"

class G4AdvancedBrepShapeRepresentationCreator: private G4GeometryCreator 
{
  public:
  
  // Constructor & destructor

    G4AdvancedBrepShapeRepresentationCreator();
    ~G4AdvancedBrepShapeRepresentationCreator();

  // Member functions

    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void* G4obj);
    const char* Name() const { return "Advanced_Brep_Shape_Representation"; }
    static G4AdvancedBrepShapeRepresentationCreator GetInstance() { return csc; }

  // Members

  private:

    static G4AdvancedBrepShapeRepresentationCreator csc;
};

#endif
