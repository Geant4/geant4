#ifndef G4ELLIPSECREATOR_HH
#define G4ELLIPSECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4EllipseCreator: private G4GeometryCreator 
{
public:
//Constructor
G4EllipseCreator();
~G4EllipseCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Ellipse";};

static G4EllipseCreator GetInstance(){return csc;};

//Members
private:

static G4EllipseCreator csc;
};
#endif
