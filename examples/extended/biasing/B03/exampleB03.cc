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
// $Id: exampleB03.cc,v 1.9 2002-11-13 08:59:22 mdressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleB03
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "G4RunManager.hh"

#include "G4CellScorerStore.hh"

// class for special output
#include "G4ScoreTable.hh"

#include "B03App.hh"


#include "Python.h"

extern "C" {
  void init_B03App();
}
int main(int argc, char **argv)
{  

  // output
  //  B03AppStruct *app = 0;
  if (argc==2 && G4String(argv[1]) == G4String("-python")) {
    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(argv[0]);
    
    /* Initialize the Python interpreter.  Required. */
    Py_Initialize();

    init_B03App();
    FILE *pyrc = fopen("python.rc","r");
    if (pyrc) {
      PyRun_SimpleFile(pyrc,"python.rc");
    }
    else {
      G4cout  << "exampleB03: python.rc not found" << G4endl;
    }
    PyRun_InteractiveLoop(stdin, "/dev/stdin");
  } 
  else {
    B03AppBase &base = B03AppBase::GetB03AppBase();
  }

  //  app = GetB03AppStruct();

  //  G4ScoreTable sp(app->fIStore);
  //  sp.Print(app->fCS_store->GetMapGeometryCellCellScorer(), &G4cout);

  return 0;
}

