#ifndef G4CONNECTEDFACESETCREATOR_HH
#define G4CONNECTEDFACESETCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4ConnectedFaceSetCreator: public G4GeometryCreator 
{
public:
//Constructor
G4ConnectedFaceSetCreator();
~G4ConnectedFaceSetCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Connected_Face_Set";};

static G4ConnectedFaceSetCreator GetInstance(){return csc;};

//Members
private:

static G4ConnectedFaceSetCreator csc;
};
#endif
