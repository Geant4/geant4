#ifndef G4CONICCREATOR_HH
#define G4CONICCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4ConicCreator: private G4GeometryCreator 
{
public:
//Constructor
G4ConicCreator();
~G4ConicCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Conic";};

static G4ConicCreator GetInstance(){return csc;};

//Members
private:

static G4ConicCreator csc;
};
#endif
