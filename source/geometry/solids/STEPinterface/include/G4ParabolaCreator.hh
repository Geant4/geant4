#ifndef G4PARABOLACREATOR_HH
#define G4PARABOLACREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4ParabolaCreator: private G4GeometryCreator 
{
public:
//Constructor
G4ParabolaCreator();
~G4ParabolaCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Parabola";};

static G4ParabolaCreator GetInstance(){return csc;};

//Members
private:

static G4ParabolaCreator csc;
};
#endif
