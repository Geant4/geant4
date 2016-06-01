#ifndef G4AXIS1PLACEMENTCREATOR_HH
#define G4AXIS1PLACEMENTCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4Axis1PlacementCreator: private G4GeometryCreator 
{
public:
//Constructor
G4Axis1PlacementCreator();
~G4Axis1PlacementCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Axis1_Placement";};

static G4Axis1PlacementCreator GetInstance(){return csc;};

//Members
private:

static G4Axis1PlacementCreator csc;
};
#endif
