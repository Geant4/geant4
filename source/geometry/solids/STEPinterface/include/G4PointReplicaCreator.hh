#ifndef G4POINTREPLICACREATOR_HH
#define G4POINTREPLICACREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4PointReplicaCreator: private G4GeometryCreator 
{
public:
//Constructor
G4PointReplicaCreator();
~G4PointReplicaCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Point_Replica";};

static G4PointReplicaCreator GetInstance(){return csc;};

//Members
private:

static G4PointReplicaCreator csc;
};
#endif
