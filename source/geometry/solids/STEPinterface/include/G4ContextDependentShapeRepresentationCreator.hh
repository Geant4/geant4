#ifndef G4CONTEXTDEPENDENTSHAPEREPRESENTATIONCREATOR_HH
#define G4CONTEXTDEPENDENTSHAPEREPRESENTATIONCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4ContextDependentShapeRepresentationCreator: private G4GeometryCreator 
{
public:
//Constructor
G4ContextDependentShapeRepresentationCreator();
~G4ContextDependentShapeRepresentationCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "Context_Dependent_Shape_Representation";};

static G4ContextDependentShapeRepresentationCreator GetInstance(){return csc;};

//Members
private:
static G4ContextDependentShapeRepresentationCreator csc;
};
#endif
