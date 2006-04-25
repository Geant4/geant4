// $Id: pyG4VisManager.cc,v 1.2 2006-04-25 08:09:46 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VisManager.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "pyG4Version.hh"
#include "G4VisManager.hh"
#if G4VERSION_NUMBER >= 800
#include "G4TrajectoryModelFactories.hh"
#endif

using namespace boost::python;

// ====================================================================
// wrappers
// ====================================================================
class PyG4VisManager : public G4VisManager {
public:
  PyG4VisManager() { SetVerboseLevel(quiet); }
  ~PyG4VisManager() { }

  void SetVerbosity(Verbosity verb) { fVerbosity= verb; }

  virtual void RegisterGraphicsSystems() { }

#if G4VERSION_NUMBER >= 800
  virtual void RegisterModelFactories() {
    RegisterModelFactory(new G4TrajectoryDrawByChargeFactory());
    RegisterModelFactory(new G4TrajectoryDrawByParticleIDFactory());
  }
#endif

};

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
    .def("SetVerbosity", &PyG4VisManager::SetVerbosity)
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

