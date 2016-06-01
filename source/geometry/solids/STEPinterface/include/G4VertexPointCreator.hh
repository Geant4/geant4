#ifndef G4VERTEXPOINTCREATOR_HH
#define G4VERTEXPOINTCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4VertexPointCreator: private G4GeometryCreator 
{
public:
//Constructor
G4VertexPointCreator();
~G4VertexPointCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Vertex_Point";};

static G4VertexPointCreator GetInstance(){return csc;};

//Members
private:

static G4VertexPointCreator csc;
};
#endif
