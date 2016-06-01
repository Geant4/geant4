#ifndef G4BOUNDEDSURFACECREATOR_HH
#define G4BOUNDEDSURFACECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"
#include "G4ControlPoints.hh"
#include "G4BSplineSurface.hh"

class G4BoundedSurfaceCreator: private G4GeometryCreator 
{
public:
//Constructor
G4BoundedSurfaceCreator();
~G4BoundedSurfaceCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "Bounded_Surface";};

static G4BoundedSurfaceCreator GetInstance(){return csc;};

//Members
private:

static G4BoundedSurfaceCreator csc;
};
#endif
