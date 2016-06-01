#ifndef G4FACEBOUNDCREATOR_HH
#define G4FACEBOUNDCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4FaceBoundCreator: public G4GeometryCreator 
{
public:
//Constructor
G4FaceBoundCreator();
~G4FaceBoundCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Face_Bound";};

static G4FaceBoundCreator GetInstance(){return csc;};

//Members
private:

static G4FaceBoundCreator csc;
};
#endif
