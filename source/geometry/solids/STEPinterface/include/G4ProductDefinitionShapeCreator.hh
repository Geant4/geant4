// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProductDefinitionShapeCreator.hh,v 1.3 2000-11-09 16:35:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ProductDefinitionShapeCreator
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
#ifndef G4PRODUCTDEFINITIONSHAPECREATOR_HH
#define G4PRODUCTDEFINITIONSHAPECREATOR_HH

#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4ProductDefinitionShapeCreator: private G4GeometryCreator 
{
  public:

  // Constructor & destructor

    G4ProductDefinitionShapeCreator();
    ~G4ProductDefinitionShapeCreator();

  // Member functions

    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void*);
    const char* Name() const { return "Product_Definition_Shape"; }
    static G4ProductDefinitionShapeCreator GetInstance() { return csc; }

  // Members

  private:

    static G4ProductDefinitionShapeCreator csc;
};

#endif
