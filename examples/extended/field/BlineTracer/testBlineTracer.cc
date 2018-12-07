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
//
// 
/// \file testBlineTracer.cc
/// \brief Test program for the G4BlineTracer class

#include "G4BlineTracer.hh"
#include "G4RunManager.hh"
#include "FTFP_BERT.hh"

// Test program which only instantiates the G4BlineTracer class.

int main()
{
  // Construct the default run manager
  // and set a physics list (otherwise a run action cannot be instatiated)
  G4RunManager * runManager = new G4RunManager;
  runManager->SetUserInitialization(new FTFP_BERT);

  // Instantiate the G4BlineTracer class
  G4BlineTracer* theBlineTool = new G4BlineTracer();
  
  // delete it
  delete theBlineTool;
}


