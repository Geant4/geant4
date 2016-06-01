#ifndef G4BSPLINESURFACEWITHKNOTSCREATOR_HH
#define G4BSPLINESURFACEWITHKNOTSCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"
#include "G4BSplineSurfaceCreator.hh"
#include "G4ControlPoints.hh"

class G4BSplineSurfaceWithKnotsCreator: public G4BSplineSurfaceCreator
{
public:
//Constructor
G4BSplineSurfaceWithKnotsCreator();
~G4BSplineSurfaceWithKnotsCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "B_Spline_Surface_With_Knots";};

static G4BSplineSurfaceWithKnotsCreator GetInstance(){return csc;};

//Members
private:

static G4BSplineSurfaceWithKnotsCreator csc;
};
#endif
