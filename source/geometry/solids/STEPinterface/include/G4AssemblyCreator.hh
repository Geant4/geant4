
#ifndef G4ASSEMBLYCREATOR_HH
#define G4ASSEMBLYCREATOR_HH

#include "G4GeometryTable.hh"
#include "G4StepFileReader.hh"
#include "G4NISTStepReader.hh"

#include "G4AdvancedBrepShapeRepresentationCreator.hh"
#include "G4ManifoldSolidBrepCreator.hh"
#include "G4ClosedShellCreator.hh"
#include "G4OpenShellCreator.hh"
#include "G4PointCreator.hh"


class G4AssemblyCreator: public G4GeometryCreator
{
public:
  // Constructor
  G4AssemblyCreator();
  G4AssemblyCreator(G4String, G4String ="NIST");

  // Destructor
  ~G4AssemblyCreator();    

  //Member functions
  void ReadStepFile();
  void CreateG4Geometry(STEPentity&);
  void CreateSTEPGeometry(void* G4obj);
  G4String Name(){return "Assembly";}
  static G4AssemblyCreator GetInstance(){return ci;}
  
  //Members
  int index; // add by L. Broglia, memory for get STEP entity
  

private:
  G4StepFileReader* StepReader;
  G4String STEPfileName;
  static G4AssemblyCreator ci;
};

#endif

