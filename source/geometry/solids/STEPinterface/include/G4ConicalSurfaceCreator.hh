#ifndef G4CONICALSURFACECREATOR_HH
#define G4CONICALSURFACECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4ConicalSurfaceCreator: private G4GeometryCreator 
{
public:
//Constructor
G4ConicalSurfaceCreator();
~G4ConicalSurfaceCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "Conical_Surface";};

static G4ConicalSurfaceCreator GetInstance(){return csc;};

//Members
private:

static G4ConicalSurfaceCreator csc;
};
#endif
