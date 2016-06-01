#ifndef G4CLOSEDSHELLCREATOR_HH
#define G4CLOSEDSHELLCREATOR_HH
#include "G4ConnectedFaceSetCreator.hh"
#include "G4GeometryTable.hh"
#include "G4BREPSolid.hh"

class G4ClosedShellCreator: public G4ConnectedFaceSetCreator 
{
public:
//Constructor
G4ClosedShellCreator();
~G4ClosedShellCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Closed_Shell";};

static G4ClosedShellCreator GetInstance(){return csc;};

//Members
private:

static G4ClosedShellCreator csc;
};
#endif
