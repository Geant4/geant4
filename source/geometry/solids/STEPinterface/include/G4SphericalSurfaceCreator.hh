#ifndef G4SPHERICALSURFACECREATOR_HH
#define G4SPHERICALSURFACECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4SphericalSurfaceCreator: private G4GeometryCreator 
{
public:
//Constructor
G4SphericalSurfaceCreator();
~G4SphericalSurfaceCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "Spherical_Surface";};

static G4SphericalSurfaceCreator GetInstance(){return csc;};

//Members
private:

static G4SphericalSurfaceCreator csc;
};
#endif
