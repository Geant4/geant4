// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ShapeRepresentationRelationshipCreator.hh,v 1.2 2000-01-21 13:45:30 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ShapeRepresentationRelationshipCreator
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
#ifndef G4SHAPEREPRESENTATIONRELATIONSHIPCREATOR_HH
#define G4SHAPEREPRESENTATIONRELATIONSHIPCREATOR_HH

#include "G4GeometryCreator.hh"

class G4ShapeRepresentationRelationshipCreator: private G4GeometryCreator 
{
  public:

  // Constructor & destructor

    G4ShapeRepresentationRelationshipCreator();
    ~G4ShapeRepresentationRelationshipCreator();

  // Member functions

    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void*);
    G4String Name() { return "Shape_Representation_Relationship"; }
    static G4ShapeRepresentationRelationshipCreator GetInstance() { return csc; }

  // Members

  private:

    static G4ShapeRepresentationRelationshipCreator csc;
};

#endif
