#ifndef G4RATIONALBSPLINESURFACECREATOR_HH
#define G4RATIONALBSPLINESURFACECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"
#include "G4BSplineSurfaceCreator.hh"

class G4RationalBSplineSurfaceCreator: public G4BSplineSurfaceCreator 
{
public:
//Constructor
G4RationalBSplineSurfaceCreator();
~G4RationalBSplineSurfaceCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "Rational_B_Spline_Surface";};

static G4RationalBSplineSurfaceCreator GetInstance(){return csc;};

//Members
private:

static G4RationalBSplineSurfaceCreator csc;
};
#endif
