#ifndef G4OPENSHELLCREATOR_HH
#define G4OPENSHELLCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4OpenShellCreator: private G4GeometryCreator 
{
public:
//Constructor
G4OpenShellCreator();
~G4OpenShellCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Open_Shell";};

static G4OpenShellCreator GetInstance(){return csc;};

//Members
private:

static G4OpenShellCreator csc;
};
#endif
