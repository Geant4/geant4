#ifndef G4AXIS2PLACEMENTCREATOR_HH
#define G4AXIS2PLACEMENTCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4Axis2PlacementCreator: private G4GeometryCreator 
{
public:
//Constructor
G4Axis2PlacementCreator();
~G4Axis2PlacementCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Axis2_Placement";};

static G4Axis2PlacementCreator GetInstance(){return csc;};

//Members
private:

static G4Axis2PlacementCreator csc;
};
#endif
