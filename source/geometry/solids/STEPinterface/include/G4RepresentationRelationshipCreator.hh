// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RepresentationRelationshipCreator.hh,v 1.4 2000-11-09 16:35:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4RepresentationRelationshipCreator
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
#ifndef G4REPRESENTATIONRELATIONSHIPCREATOR_HH
#define G4REPRESENTATIONRELATIONSHIPCREATOR_HH

#include "G4GeometryCreator.hh"

class G4RepresentationRelationshipCreator: private G4GeometryCreator 
{
  public:

  // Constructor & destructor

    G4RepresentationRelationshipCreator();
    ~G4RepresentationRelationshipCreator();

  // Member functions

    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void*);
    const char* Name() const { return "Representation_Relationship"; }
    static G4RepresentationRelationshipCreator GetInstance() { return csc; }

  // Members

  private:

    static G4RepresentationRelationshipCreator csc;
};

#endif
