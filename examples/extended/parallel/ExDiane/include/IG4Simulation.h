#ifndef IG4SIMULATION_H
#define IG4SIMULATION_H

/* ==========================================
   DIANE - Distributed Analysis Environment 
   Copyright (C) Jakub T. Moscicki, 2000-2003
   ------------------------------------------
   See $DIANE_TOP/LICENSE for details.
   ========================================== 
*/

# include <string>

namespace DIANE
{

class IG4Simulation
{

public:
  
  virtual void setSeed(int seed) = 0;  
  virtual bool initialize(int argc, char** argv) = 0;
  virtual void executeMacro(std::string macroFileName) = 0;
  virtual std::string getOutputFilename() = 0;

  virtual ~IG4Simulation() {}
};
}

extern "C" 
DIANE::IG4Simulation* createG4Simulation(int);

#endif

