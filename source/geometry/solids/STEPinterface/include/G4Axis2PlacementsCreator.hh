#ifndef G4AXIS2PLACEMENTSCREATOR_HH
#define G4AXIS2PLACEMENTSCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4Axis2PlacementsCreator: private G4GeometryCreator 
{
public:
//Constructor
G4Axis2PlacementsCreator();
~G4Axis2PlacementsCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Axis2_Placements";};

static G4Axis2PlacementsCreator GetInstance(){return csc;};

//Members
private:

static G4Axis2PlacementsCreator csc;
};
#endif
