#ifndef G4POINTONSURFACECREATOR_HH
#define G4POINTONSURFACECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4PointOnSurfaceCreator: private G4GeometryCreator 
{
public:
//Constructor
G4PointOnSurfaceCreator();
~G4PointOnSurfaceCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Point_On_Surface";};

static G4PointOnSurfaceCreator GetInstance(){return csc;};

//Members
private:

static G4PointOnSurfaceCreator csc;
};
#endif
