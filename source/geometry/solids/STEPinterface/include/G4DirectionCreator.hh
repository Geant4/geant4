#ifndef G4DIRECTIONCREATOR_HH
#define G4DIRECTIONCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4DirectionCreator: private G4GeometryCreator 
{
public:
//Constructor
G4DirectionCreator();
~G4DirectionCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Direction";};

static G4DirectionCreator GetInstance(){return csc;};

//Members
private:

static G4DirectionCreator csc;
};
#endif
