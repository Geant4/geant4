#ifndef G4SHAPEREPRESENTATIONCREATOR_HH
#define G4SHAPEREPRESENTATIONCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"
#include "G4PlacementVector.hh"
class G4ShapeRepresentationCreator: private G4GeometryCreator 
{
public:
//Constructor
G4ShapeRepresentationCreator();
~G4ShapeRepresentationCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "Shape_Representation";};

static G4ShapeRepresentationCreator GetInstance(){return csc;};

//Members
private:

static G4ShapeRepresentationCreator csc;
};
#endif
