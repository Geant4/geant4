#include "globals.hh"

#include "Randomize.hh"
#include "IG4Simulation.h"
#include "G4RunManager.hh"
class BrachySimulation : virtual public DIANE::IG4Simulation 

{
public:
  BrachySimulation(int);
  ~BrachySimulation();
  
   void setSeed(int seed);  
   bool initialize(int argc, char** argv);
   void executeMacro(std::string macroFileName);
   std::string getOutputFilename();

private: 
  int seed;
};



