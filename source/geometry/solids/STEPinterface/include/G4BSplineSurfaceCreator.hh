#ifndef G4BSPLINESURFACECREATOR_HH
#define G4BSPLINESURFACECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4BSplineSurfaceCreator: public G4GeometryCreator 
{
public:
//Constructor
G4BSplineSurfaceCreator();
~G4BSplineSurfaceCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "B_Spline_Surface";};

static G4BSplineSurfaceCreator GetInstance(){return csc;};

//Members
private:

static G4BSplineSurfaceCreator csc;
};
#endif
