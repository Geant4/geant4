#ifndef G4VECTORCREATOR_HH
#define G4VECTORCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4VectorCreator: private G4GeometryCreator 
{
public:
//Constructor
G4VectorCreator();
~G4VectorCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Vector";};

static G4VectorCreator GetInstance(){return csc;};

//Members
private:

static G4VectorCreator csc;
};
#endif
