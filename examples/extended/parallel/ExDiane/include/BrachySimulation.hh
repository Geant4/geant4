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
//
// $Id: BrachySimulation.hh
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// --------------------------------------------------------------
//                 GEANT 4 - Brachytherapy example
// --------------------------------------------------------------
//
// Code developed by: S.Guatelli, K. Moscicki
//
//    *******************************
//    *                             *
//    *    BrachySimulation.cc      *
//    *                             *
//    *******************************
//

#include "globals.hh"

#include "Randomize.hh"
#include "IG4Simulation.h"
#include "G4RunManager.hh"
class BrachySimulation : virtual public DIANE::IG4Simulation 
{
public:
  BrachySimulation(G4int);
  ~BrachySimulation();
  
  void setSeed(G4int seed);  
  G4bool initialize(int argc, char** argv);
  void executeMacro(std::string macroFileName);
  std::string getOutputFilename();
  void finish();

private: 
  G4int seed; 
  G4RunManager* pRunManager;
};



