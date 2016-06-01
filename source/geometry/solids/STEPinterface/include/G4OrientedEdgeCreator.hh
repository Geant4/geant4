#ifndef G4ORIENTEDEDGECREATOR_HH
#define G4ORIENTEDEDGECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4OrientedEdgeCreator: private G4GeometryCreator 
{
public:
//Constructor
G4OrientedEdgeCreator();
~G4OrientedEdgeCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Oriented_Edge";};

static G4OrientedEdgeCreator GetInstance(){return csc;};

//Members
private:

static G4OrientedEdgeCreator csc;
};
#endif
