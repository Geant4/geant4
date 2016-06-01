#ifndef G4ADVANCEDBREPSHAPEREPRESENTATIONCREATOR_HH
#define G4ADVANCEDBREPSHAPEREPRESENTATIONCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4AdvancedBrepShapeRepresentationCreator: private G4GeometryCreator 
{
public:
//Constructor
G4AdvancedBrepShapeRepresentationCreator();
~G4AdvancedBrepShapeRepresentationCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Advanced_Brep_Shape_Representation";};

static G4AdvancedBrepShapeRepresentationCreator GetInstance(){return csc;};

//Members
private:

static G4AdvancedBrepShapeRepresentationCreator csc;
};
#endif
