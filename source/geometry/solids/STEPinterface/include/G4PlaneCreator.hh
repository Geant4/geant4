#ifndef G4PLANECREATOR_HH
#define G4PLANECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4PlaneCreator: private G4GeometryCreator 
{
public:
//Constructor
G4PlaneCreator();
~G4PlaneCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Plane";};

static G4PlaneCreator GetInstance(){return csc;};

//Members
private:

static G4PlaneCreator csc;
};
#endif
