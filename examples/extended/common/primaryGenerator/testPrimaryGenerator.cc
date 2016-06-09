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
// $Id$
// 
/// \file testPrimaryGenerator.cc
/// \brief Test program for the primaryGenerator common classes

#include "ExG4PrimaryGeneratorAction01.hh"
#include "ExG4PrimaryGeneratorAction02.hh"

#include "G4RunManager.hh"
#include "QGSP_BERT.hh"

// test program which only instantiates classes defined in 
// examples/common/primaryGenerator 

int main()
{
  // First construct necessary classes
  //
  G4RunManager * runManager = new G4RunManager;
  G4VModularPhysicsList* physicsList = new QGSP_BERT;
  runManager->SetUserInitialization(physicsList);

  // Instantiate all primary generator actions classes
  ExG4PrimaryGeneratorAction01* primaryGeneratorAction01
    = new ExG4PrimaryGeneratorAction01();
  ExG4PrimaryGeneratorAction02* primaryGeneratorAction02
    = new ExG4PrimaryGeneratorAction02();
  
  // delete all
  delete primaryGeneratorAction01;
  delete primaryGeneratorAction02;

  return 0;
}


