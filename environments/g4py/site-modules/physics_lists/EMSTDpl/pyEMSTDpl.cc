// $Id: pyEMSTDpl.cc,v 1.2 2006-04-25 10:31:40 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyEMSTDpl.cc
//
//   [Qgeom]
//   a site-module of Geant4Py
//
//   Electron/Gamma EM-standard physics list
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4RunManager.hh"
#include "PhysicsListEMstd.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyEMSTDpl {

PhysicsListEMstd* Construct()
{
  PhysicsListEMstd* pl= new PhysicsListEMstd;
  
  G4RunManager* runMgr= G4RunManager::GetRunManager();
  runMgr-> SetUserInitialization(pl);

  return pl;
}

};

using namespace pyEMSTDpl;

// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(EMSTDpl) {

  class_<PhysicsListEMstd, PhysicsListEMstd*, bases<G4VUserPhysicsList> >
    ("PhysicsListEMstd", "Electron/Gamma EM-standard physics list")
    ;

  // ---
  def("Construct", Construct, 
      return_value_policy<reference_existing_object>());

}
