// $Id: pyExN03pl.cc,v 1.1 2006-04-25 10:32:41 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyExN03pl.cc
//
//   [ExN03pl]
//   a site-module of Geant4Py
//
//   a physics list presnted in ExN03 of Geant4 example
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4RunManager.hh"
#include "ExN03PhysicsList.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyExN03pl {

ExN03PhysicsList* Construct()
{
  ExN03PhysicsList* pl= new ExN03PhysicsList;

  G4RunManager* runMgr= G4RunManager::GetRunManager();
  runMgr-> SetUserInitialization(pl);

  return pl;
}

};

using namespace pyExN03pl;

// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(ExN03pl) {

  class_<ExN03PhysicsList, ExN03PhysicsList*, bases<G4VUserPhysicsList> >
    ("ExN03PhysicsList", "ExN03 physics list")
    ;

  // ---
  def("Construct", Construct,
      return_value_policy<reference_existing_object>());
  
}
