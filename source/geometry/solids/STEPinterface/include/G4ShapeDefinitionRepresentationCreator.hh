// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ShapeDefinitionRepresentationCreator.hh,v 1.2 2000-01-21 13:45:30 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ShapeDefinitionRepresentationCreator
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
#ifndef G4SHAPEDEFINITIONREPRESENTATIONCREATOR_HH
#define G4SHAPEDEFINITIONREPRESENTATIONCREATOR_HH

#include "G4GeometryCreator.hh"

class G4ShapeDefinitionRepresentationCreator: private G4GeometryCreator 
{
  public:

  // Constructor & destructor

    G4ShapeDefinitionRepresentationCreator();
    ~G4ShapeDefinitionRepresentationCreator();

  // Member functions

    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void*);
    G4String Name() { return "Shape_Definition_Representation"; }
    static G4ShapeDefinitionRepresentationCreator GetInstance() { return csc; }

  // Members

  private:

    static G4ShapeDefinitionRepresentationCreator csc;
};

#endif
