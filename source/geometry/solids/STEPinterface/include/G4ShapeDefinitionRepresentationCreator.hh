#ifndef G4SHAPEDEFINITIONREPRESENTATIONCREATOR_HH
#define G4SHAPEDEFINITIONREPRESENTATIONCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4ShapeDefinitionRepresentationCreator: private G4GeometryCreator 
{
public:
//Constructor
G4ShapeDefinitionRepresentationCreator();
~G4ShapeDefinitionRepresentationCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "Shape_Definition_Representation";};

static G4ShapeDefinitionRepresentationCreator GetInstance(){return csc;};

//Members
private:

static G4ShapeDefinitionRepresentationCreator csc;
};
#endif
