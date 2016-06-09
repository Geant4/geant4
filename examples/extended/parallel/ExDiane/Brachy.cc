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
// $Id: Brachy.cc
// GEANT4 tag $Name: geant4-07-01 $
//
// --------------------------------------------------------------
//                 GEANT 4 - ExDiane example
// --------------------------------------------------------------
//
// Code developed by: S.Guatelli

#include "BrachySimulation.hh"
int main(int argc, char** argv){

  BrachySimulation * simulation = new BrachySimulation(0);  

  simulation -> initialize(argc,argv);

  G4String mn;

  if(argc>1)
    mn = argv[1];

  if(mn=="")
    {
      std::cerr << "error: please give a macro file name" << std::endl;
      return -1;
    }
  
  simulation -> executeMacro(mn);
 
  simulation -> finish();
 
  delete simulation;
  G4cout << "The job ends! " << G4endl;
  return 0;
}
