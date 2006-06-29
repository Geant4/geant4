//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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

