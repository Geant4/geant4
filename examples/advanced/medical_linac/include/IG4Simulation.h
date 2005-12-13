//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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

