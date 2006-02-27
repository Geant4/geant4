// $Id: pyEMSTDpl.cc,v 1.1 2006-02-27 09:50:13 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyEMSTDpl.cc
//
//   [Qgeom]
//   a site-module of Geant4Py
//
//   EM-std physics list
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4RunManager.hh"
#include "QPhysicsList.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyEMSTDpl {

void Construct()
{
  G4RunManager* runMgr= G4RunManager::GetRunManager();
  runMgr-> SetUserInitialization(new QPhysicsList);
}

};

using namespace pyEMSTDpl;

// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(EMSTDpl) {

  class_<QPhysicsList, QPhysicsList*,
    bases<G4VUserPhysicsList> >
    ("Qem_stdPhysicsList", "EM-std physics list")
    ;

  // ---
  def("Construct",  Construct);
}

