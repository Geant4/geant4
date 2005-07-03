#include "globals.hh"

#include "Randomize.hh"
#include "IG4Simulation.h"
#include "G4RunManager.hh"
class MedLinacSimulation : virtual public DIANE::IG4Simulation 

{

public:
  
  MedLinacSimulation(int);
  
  ~MedLinacSimulation();
  G4RunManager* pRunManager; 
   void setSeed(int seed);  
   bool initialize(int argc, char** argv);
   void executeMacro(std::string macroFileName);
   std::string getOutputFilename();

private: 

  int seed;

};



