#ifndef G4BSPLINECURVEWITHKNOTSCREATOR_HH
#define G4BSPLINECURVEWITHKNOTSCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4BSplineCurveWithKnotsCreator: private G4GeometryCreator 
{
public:
//Constructor
G4BSplineCurveWithKnotsCreator();
~G4BSplineCurveWithKnotsCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "B_Spline_Curve_With_Knots";};

static G4BSplineCurveWithKnotsCreator GetInstance(){return csc;};

//Members
private:

static G4BSplineCurveWithKnotsCreator csc;
};
#endif
