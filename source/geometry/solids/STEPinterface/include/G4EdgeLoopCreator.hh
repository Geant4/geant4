#ifndef G4EDGELOOPCREATOR_HH
#define G4EDGELOOPCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4EdgeLoopCreator: private G4GeometryCreator 
{
public:
//Constructor
G4EdgeLoopCreator();
~G4EdgeLoopCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Edge_Loop";};

static G4EdgeLoopCreator GetInstance(){return csc;};

//Members
private:

static G4EdgeLoopCreator csc;
};
#endif
