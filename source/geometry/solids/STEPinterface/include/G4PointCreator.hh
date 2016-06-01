#ifndef G4POINTCREATOR_HH
#define G4POINTCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4PointCreator: private G4GeometryCreator 
{
public:
//Constructor
G4PointCreator();
~G4PointCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Point";};

static G4PointCreator GetInstance(){return csc;};

//Members
private:

static G4PointCreator csc;
};
#endif
