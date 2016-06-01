#ifndef G4CARTESIANPOINTCREATOR_HH
#define G4CARTESIANPOINTCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4CartesianPointCreator: private G4GeometryCreator 
{
public:
//Constructor
G4CartesianPointCreator();
~G4CartesianPointCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Cartesian_Point";};

static G4CartesianPointCreator GetInstance(){return csc;};

//Members
private:

static G4CartesianPointCreator csc;
};
#endif
