#ifndef G4CIRCLECREATOR_HH
#define G4CIRCLECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4CircleCreator: private G4GeometryCreator 
{
public:
//Constructor
G4CircleCreator();
~G4CircleCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Circle";};

static G4CircleCreator GetInstance(){return csc;};

//Members
private:

static G4CircleCreator csc;
};
#endif
