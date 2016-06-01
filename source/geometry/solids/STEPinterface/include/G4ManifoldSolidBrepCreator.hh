#ifndef G4MANIFOLDSOLIDBREPCREATOR_HH
#define G4MANIFOLDSOLIDBREPCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"
#include "G4BREPSolid.hh"

class G4ManifoldSolidBrepCreator: private G4GeometryCreator 
{
public:
//Constructor
G4ManifoldSolidBrepCreator();
~G4ManifoldSolidBrepCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Manifold_Solid_Brep";};

static G4ManifoldSolidBrepCreator GetInstance(){return csc;};

//Members
private:

static G4ManifoldSolidBrepCreator csc;
};
#endif
