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
// $Id: pymodG4run.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
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
