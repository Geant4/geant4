#ifndef G4CURVECREATOR_HH
#define G4CURVECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4CurveCreator: private G4GeometryCreator 
{
public:
//Constructor
G4CurveCreator();
~G4CurveCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Curve";};

static G4CurveCreator GetInstance(){return csc;};

//Members
private:

static G4CurveCreator csc;
};
#endif
