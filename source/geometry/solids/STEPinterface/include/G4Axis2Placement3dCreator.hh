#ifndef G4AXIS2PLACEMENT3DCREATOR_HH
#define G4AXIS2PLACEMENT3DCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4Axis2Placement3dCreator: private G4GeometryCreator 
{
public:
//Constructor
G4Axis2Placement3dCreator();
~G4Axis2Placement3dCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "Axis2_Placement_3d";};

static G4Axis2Placement3dCreator GetInstance(){return csc;};

//Members
private:

static G4Axis2Placement3dCreator csc;
};
#endif
