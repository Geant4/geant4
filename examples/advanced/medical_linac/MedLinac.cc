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
// $Id: MedLinac.cc,v 1.6 2005-07-03 23:27:36 mpiergen Exp $
//
// --------------------------------------------------------------
//      GEANT 4 -  medical_linac
//
// Code developed by: M. Piergentili

#include "MedLinacSimulation.hh"

int main(int argc, char** argv){

  MedLinacSimulation * simulation = new MedLinacSimulation(0);  

  simulation->initialize(argc,argv);

  std::string mn;

  if(argc>1)
    mn = argv[1];

  if(mn=="")
    {
      std::cerr << "error: please give a macro file name" << std::endl;
      return -1;
    }
  
  simulation->executeMacro(mn);
  
  delete simulation;

  return 0;
}
