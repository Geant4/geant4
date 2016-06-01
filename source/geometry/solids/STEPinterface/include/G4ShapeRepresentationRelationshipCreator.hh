#ifndef G4SHAPEREPRESENTATIONRELATIONSHIPCREATOR_HH
#define G4SHAPEREPRESENTATIONRELATIONSHIPCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"
#include "G4PlacementVector.hh"

class G4ShapeRepresentationRelationshipCreator: private G4GeometryCreator 
{
public:
//Constructor
G4ShapeRepresentationRelationshipCreator();
~G4ShapeRepresentationRelationshipCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "Shape_Representation_Relationship";};

static G4ShapeRepresentationRelationshipCreator GetInstance(){return csc;};

//Members
private:

static G4ShapeRepresentationRelationshipCreator csc;
};
#endif
