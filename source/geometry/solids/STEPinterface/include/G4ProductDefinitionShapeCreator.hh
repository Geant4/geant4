#ifndef G4PRODUCTDEFINITIONSHAPECREATOR_HH
#define G4PRODUCTDEFINITIONSHAPECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4ProductDefinitionShapeCreator: private G4GeometryCreator 
{
public:
//Constructor
G4ProductDefinitionShapeCreator();
~G4ProductDefinitionShapeCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "Product_Definition_Shape";};

static G4ProductDefinitionShapeCreator GetInstance(){return csc;};

//Members
private:

static G4ProductDefinitionShapeCreator csc;
};
#endif
