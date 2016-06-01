#ifndef G4HYPERBOLACREATOR_HH
#define G4HYPERBOLACREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4HyperbolaCreator: private G4GeometryCreator 
{
public:
//Constructor
G4HyperbolaCreator();
~G4HyperbolaCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Hyperbola";};

static G4HyperbolaCreator GetInstance(){return csc;};

//Members
private:

static G4HyperbolaCreator csc;
};
#endif
