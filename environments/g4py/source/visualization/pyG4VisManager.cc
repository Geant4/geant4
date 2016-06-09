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
// $Id: pyG4VisManager.cc,v 1.7 2008-03-13 07:32:18 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VisManager.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VisManager.hh"
#include "G4TrajectoryModelFactories.hh"

using namespace boost::python;

// ====================================================================
// wrappers
// ====================================================================
class PyG4VisManager : public G4VisManager {
public:
  PyG4VisManager() { SetVerboseLevel(quiet); }
  ~PyG4VisManager() { }

  virtual void RegisterGraphicsSystems() { }

  virtual void RegisterModelFactories() {
    RegisterModelFactory(new G4TrajectoryDrawByChargeFactory());
    RegisterModelFactory(new G4TrajectoryDrawByParticleIDFactory());
  }

};

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4VisManager {

void (PyG4VisManager::*f1_SetVerboseLevel)(G4int)
  = &PyG4VisManager::SetVerboseLevel;
void (PyG4VisManager::*f2_SetVerboseLevel)(const G4String&)
  = &PyG4VisManager::SetVerboseLevel;
  void (PyG4VisManager::*f3_SetVerboseLevel)(G4VisManager::Verbosity)
  = &PyG4VisManager::SetVerboseLevel;

}

using namespace pyG4VisManager;

// ====================================================================
// module definition
// ====================================================================
void export_G4VisManager()
{
  scope in_PyG4VisManager =
    class_<PyG4VisManager, boost::noncopyable>
    ("G4VisManager", "visualization manager")
    // ---
    .def("GetConcreteInstance", &PyG4VisManager::GetConcreteInstance,
         "Get an instance of G4VisManager",
         return_value_policy<reference_existing_object>())
    .staticmethod("GetConcreteInstance")
    // ---
    .def("SetVerboseLevel", f1_SetVerboseLevel)
    .def("SetVerboseLevel", f2_SetVerboseLevel)
    .def("SetVerboseLevel", f3_SetVerboseLevel)
    .def("GetVerbosity", &PyG4VisManager::GetVerbosity)
    .def("Initialize", &PyG4VisManager::Initialize)
    .def("RegisterGraphicsSystem", &PyG4VisManager::RegisterGraphicsSystem)
    ;

  // enum LineStyle
  enum_<G4VisManager::Verbosity>("Verbosity")
    .value("quiet",           G4VisManager::quiet)
    .value("startuo",         G4VisManager::startup)
    .value("errors",          G4VisManager::errors)
    .value("warnings",        G4VisManager::warnings)
    .value("confirmations",   G4VisManager::confirmations)
    .value("parameters",      G4VisManager::parameters)
    .value("all",             G4VisManager::all)
    ;
}

