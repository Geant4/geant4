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
// $Id: pyG4RunManagerKernel.cc 86749 2014-11-17 15:03:05Z gcosmo $
// ====================================================================
//   pyG4RunManagerKernel.cc
//
//                                         2006 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4RunManagerKernel.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4RunManagerKernel {

// RunInitialization()
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_RunInitialization,
                                       RunInitialization, 0, 1)
    
}

using namespace pyG4RunManagerKernel;

// ====================================================================
// module definition
// ====================================================================
void export_G4RunManagerKernel()
{
  class_<G4RunManagerKernel>("G4RunManagerKernel", "run manager kernel")
    .def("GetRunManagerKernel", &G4RunManagerKernel::GetRunManagerKernel,
         "Get an instance of G4RunManagerKernel",
         return_value_policy<reference_existing_object>())
    .staticmethod("GetRunManagerKernel")
    // ---
    //.def("DefineWorldVolume", &G4RunManagerKernel::DefineWorldVolume)
    //.def("SetPhysics", &G4RunManagerKernel::SetPhysics)
    //.def("InitializePhysics", &G4RunManagerKernel::InitializePhysics)
    .def("RunInitialization",  &G4RunManagerKernel::RunInitialization,
         f_RunInitialization())
    //.def("RunTermination", &G4RunManagerKernel::RunTermination)
    //.def("UpdateRegion", &G4RunManagerKernel::UpdateRegion)
    //.def("DumpRegion", &G4RunManagerKernel::DumpRegion)
    //.def("DumpRegion", &G4RunManagerKernel::DumpRegion)
    //.def("GeometryHasBeenModified",
    //&G4RunManagerKernel::GeometryHasBeenModified)
    //.def("PhysicsHasBeenModified",
    //&G4RunManagerKernel::PhysicsHasBeenModified)
    //.def("GetEventManager", &G4RunManagerKernel::GetEventManager,
    //...)
    //.def("GetStackManager", &G4RunManagerKernel::GetStackManager,
    //...)
    //.def("GetTrackingManager", &G4RunManagerKernel::GetTrackingManager,
    //...)
    //.def("SetPrimaryTransformer", &G4RunManagerKernel::SetPrimaryTransformer)
    //.def("GetPrimaryTransformer", &G4RunManagerKernel::GetPrimaryTransformer,
    //...)
    //.def("GetVersionString", &G4RunManagerKernel::GetVersionString)
    //.def("SetVerboseLevel", &G4RunManagerKernel::SetVerboseLevel)
    //.def("SetGeometryToBeOptimized",
    //&G4RunManagerKernel::SetGeometryToBeOptimized)
    ;
}
