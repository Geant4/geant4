#ifndef G4GEOMETRICREPRESENTATIONCONTEXTCREATOR_HH
#define G4GEOMETRICREPRESENTATIONCONTEXTCREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4GeometricRepresentationContextCreator: private G4GeometryCreator 
{
public:
//Constructor
G4GeometricRepresentationContextCreator();
~G4GeometricRepresentationContextCreator();

//Member functions

void CreateG4Geometry(STEPentity&);

void CreateSTEPGeometry(void*);

G4String Name(){return "Geometric_Representation_Context";};

static G4GeometricRepresentationContextCreator GetInstance(){return csc;};

//Members
private:

static G4GeometricRepresentationContextCreator csc;
};
#endif
