#ifndef G4TOROIDALSURFACECREATOR_HH
#define G4TOROIDALSURFACECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4ToroidalSurfaceCreator: private G4GeometryCreator 
{
public:
//Constructor
G4ToroidalSurfaceCreator();
~G4ToroidalSurfaceCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Toroidal_Surface";};

static G4ToroidalSurfaceCreator GetInstance(){return csc;};

//Members
private:

static G4ToroidalSurfaceCreator csc;
};
#endif
