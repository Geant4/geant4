#ifndef G4LINECREATOR_HH
#define G4LINECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"
#include "G4Line.hh"

class G4LineCreator: private G4GeometryCreator 
{
public:
//Constructor
G4LineCreator();
~G4LineCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Line";};

static G4LineCreator GetInstance(){return csc;};

//Members
private:

static G4LineCreator csc;
};
#endif
