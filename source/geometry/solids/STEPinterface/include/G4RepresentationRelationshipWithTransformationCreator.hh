// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RepresentationRelationshipWithTransformationCreator.hh,v 1.3 2000-11-09 16:35:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4RepresentationRelationshipWithTransformationCreator
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
#ifndef G4REPRESENTATIONRELATIONSHIPWITHTRANSFORMATIONCREATOR_HH
#define G4REPRESENTATIONRELATIONSHIPWITHTRANSFORMATIONCREATOR_HH

#include "G4GeometryCreator.hh"

class G4RepresentationRelationshipWithTransformationCreator
      : private G4GeometryCreator 
{
  public:
  
  // Constructor & destructor

    G4RepresentationRelationshipWithTransformationCreator();
    ~G4RepresentationRelationshipWithTransformationCreator();

  // Member functions

    void CreateG4Geometry(STEPentity&);
    void CreateSTEPGeometry(void*);
    const char* Name() const
      { return "Representation_Relationship_With_Transformation"; }
    static G4RepresentationRelationshipWithTransformationCreator GetInstance()
      { return csc; }

  // Members

  private:

    static G4RepresentationRelationshipWithTransformationCreator csc;
};

#endif
