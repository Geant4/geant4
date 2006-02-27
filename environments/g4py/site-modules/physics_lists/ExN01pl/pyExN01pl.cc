// $Id: pyExN01pl.cc,v 1.1 2006-02-27 09:50:13 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyExN01pl.cc
//
//   [ExN01pl]
//   a site-module of Geant4Py
//
//   minimal physics list, which is presnted in ExN01 of Geant4 example
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4RunManager.hh"
#include "ExN01PhysicsList.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyExN01pl {

void Construct()
{
  G4RunManager* runMgr= G4RunManager::GetRunManager();
  runMgr-> SetUserInitialization(new ExN01PhysicsList);
}

};

using namespace pyExN01pl;

// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(ExN01pl) {

  class_<ExN01PhysicsList, ExN01PhysicsList*,
    bases<G4VUserPhysicsList> >
    ("ExN01PhysicsList", "ExN01 physics list")
    ;

  // ---
  def("Construct",  Construct);
}

