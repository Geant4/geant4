#ifndef G4ADVANCEDFACECREATOR_HH
#define G4ADVANCEDFACECREATOR_HH
#include "G4GeometryCreator.hh"
#include "G4GeometryTable.hh"

class G4AdvancedFaceCreator;  // Forward declararation -> Simone

class G4AdvancedFaceCreator: private G4GeometryCreator 
{
public:
//Constructor
G4AdvancedFaceCreator();
~G4AdvancedFaceCreator();

//Member functions
void CreateG4Geometry(STEPentity&);
void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Advanced_Face";};

static G4AdvancedFaceCreator GetInstance(){return csc;};

//Members
private:

static G4AdvancedFaceCreator csc;
};
#endif
