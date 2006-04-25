// $Id: pyExN01pl.cc,v 1.2 2006-04-25 10:31:40 kmura Exp $
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

ExN01PhysicsList* Construct()
{
  ExN01PhysicsList* pl= new ExN01PhysicsList;

  G4RunManager* runMgr= G4RunManager::GetRunManager();
  runMgr-> SetUserInitialization(pl);

  return pl;
}

};

using namespace pyExN01pl;

// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(ExN01pl) {

  class_<ExN01PhysicsList, ExN01PhysicsList*, bases<G4VUserPhysicsList> >
    ("ExN01PhysicsList", "ExN01 physics list")
    ;

  // ---
  def("Construct", Construct,
      return_value_policy<reference_existing_object>());
  
}
