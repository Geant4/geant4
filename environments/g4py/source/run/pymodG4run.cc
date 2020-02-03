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
// ====================================================================
//   pymodG4run.cc [Geant4Py module]
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================

void export_G4RunManager();
void export_G4RunManagerKernel();
void export_G4Run();
void export_G4UserRunAction();
void export_G4VUserPrimaryGeneratorAction();
void export_G4VUserDetectorConstruction();
void export_G4VUserPhysicsList();
void export_G4VModularPhysicsList();
void export_G4VPhysicsConstructor();

BOOST_PYTHON_MODULE(G4run) 
{
  export_G4RunManager();
  export_G4RunManagerKernel();
  export_G4Run();
  export_G4UserRunAction();
  export_G4VUserPrimaryGeneratorAction();
  export_G4VUserDetectorConstruction();
  export_G4VUserPhysicsList();
  export_G4VModularPhysicsList();
  export_G4VPhysicsConstructor();
}
