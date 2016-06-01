#ifndef G4POINTONCURVECREATOR_HH
#define G4POINTONCURVECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4PointOnCurveCreator: private G4GeometryCreator 
{
public:
//Constructor
G4PointOnCurveCreator();
~G4PointOnCurveCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Point_On_Curve";};

static G4PointOnCurveCreator GetInstance(){return csc;};

//Members
private:

static G4PointOnCurveCreator csc;
};
#endif
