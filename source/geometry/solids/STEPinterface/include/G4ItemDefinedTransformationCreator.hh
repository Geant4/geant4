#ifndef G4ITEMDEFINEDTRANSFORMATIONCREATOR_HH
#define G4ITEMDEFINEDTRANSFORMATIONCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"
#include "G4PlacementVector.hh"

class G4ItemDefinedTransformationCreator: private G4GeometryCreator 
{
public:
//Constructor
G4ItemDefinedTransformationCreator();
~G4ItemDefinedTransformationCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "Item_Defined_Transformation";};

static G4ItemDefinedTransformationCreator GetInstance(){return csc;};

//Members
private:

static G4ItemDefinedTransformationCreator csc;
};
#endif
