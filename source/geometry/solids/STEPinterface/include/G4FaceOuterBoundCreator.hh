#ifndef G4FACEOUTERBOUNDCREATOR_HH
#define G4FACEOUTERBOUNDCREATOR_HH
#include "G4GeometryCreator.hh"
//#include "G4GeometryTable.hh"
#include "G4FaceBoundCreator.hh"

class G4FaceOuterBoundCreator: public G4FaceBoundCreator
{
public:
//Constructor
G4FaceOuterBoundCreator();
~G4FaceOuterBoundCreator();

//Member functions

void CreateSTEPGeometry(void* G4obj);

G4String Name(){return "Face_Outer_Bound";};

static G4FaceOuterBoundCreator GetInstance(){return csc;};

//Members
private:

static G4FaceOuterBoundCreator csc;
};
#endif


