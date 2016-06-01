#ifndef G4AXIS2PLACEMENT2DCREATOR_HH
#define G4AXIS2PLACEMENT2DCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4Axis2Placement2dCreator: private G4GeometryCreator 
{
public:
//Constructor
G4Axis2Placement2dCreator();
~G4Axis2Placement2dCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Axis2_Placement_2d";};

static G4Axis2Placement2dCreator GetInstance(){return csc;};

//Members
private:

static G4Axis2Placement2dCreator csc;
};
#endif
