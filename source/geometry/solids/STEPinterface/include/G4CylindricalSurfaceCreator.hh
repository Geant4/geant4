#ifndef G4CYLINDRICALSURFACECREATOR_HH
#define G4CYLINDRICALSURFACECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4CylindricalSurfaceCreator: private G4GeometryCreator 
{
public:
//Constructor
G4CylindricalSurfaceCreator();
~G4CylindricalSurfaceCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Cylindrical_Surface";};

static G4CylindricalSurfaceCreator GetInstance(){return csc;};

//Members
private:

static G4CylindricalSurfaceCreator csc;
};
#endif
