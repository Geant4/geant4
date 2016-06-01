#ifndef G4REPRESENTATIONRELATIONSHIPCREATOR_HH
#define G4REPRESENTATIONRELATIONSHIPCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4RepresentationRelationshipCreator: private G4GeometryCreator 
{
public:
//Constructor
G4RepresentationRelationshipCreator();
~G4RepresentationRelationshipCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "Representation_Relationship";};

static G4RepresentationRelationshipCreator GetInstance(){return csc;};

//Members
private:
static G4int placeCount;
static G4RepresentationRelationshipCreator csc;
};
#endif
