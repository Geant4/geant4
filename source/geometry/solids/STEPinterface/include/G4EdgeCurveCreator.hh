#ifndef G4EDGECURVECREATOR_HH
#define G4EDGECURVECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4EdgeCurveCreator: private G4GeometryCreator 
{
public:
//Constructor
G4EdgeCurveCreator();
~G4EdgeCurveCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Edge_Curve";};

static G4EdgeCurveCreator GetInstance(){return csc;};

//Members
private:

static G4EdgeCurveCreator csc;
};
#endif
