#ifndef G4REPRESENTATIONRELATIONSHIPWITHTRANSFORMATIONCREATOR_HH
#define G4REPRESENTATIONRELATIONSHIPWITHTRANSFORMATIONCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4RepresentationRelationshipWithTransformationCreator: private G4GeometryCreator 
{
public:
//Constructor
G4RepresentationRelationshipWithTransformationCreator();
~G4RepresentationRelationshipWithTransformationCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "Representation_Relationship_With_Transformation";};

static G4RepresentationRelationshipWithTransformationCreator GetInstance(){return csc;};

//Members
private:

static G4RepresentationRelationshipWithTransformationCreator csc;
};
#endif
