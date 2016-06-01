#ifndef G4FACESURFACECREATOR_HH
#define G4FACESURFACECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4FaceSurfaceCreator: public G4GeometryCreator 
{
public:
//Constructor
G4FaceSurfaceCreator();
~G4FaceSurfaceCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Face_Surface";};

static G4FaceSurfaceCreator GetInstance(){return csc;};

//Members
private:

static G4FaceSurfaceCreator csc;
};
#endif
