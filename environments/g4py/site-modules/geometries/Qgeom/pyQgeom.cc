// $Id: pyQgeom.cc,v 1.2 2006-04-25 10:25:12 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyQgeom.cc
//
//   [Qgeom]
//   a site-module of Geant4Py
//
//   A sample geometry
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4RunManager.hh"
#include "QDetectorConstruction.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyQgeom {

void Construct()
{
  G4RunManager* runMgr= G4RunManager::GetRunManager();
  runMgr-> SetUserInitialization(new QDetectorConstruction);
}

};

using namespace pyQgeom;

// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(Qgeom) {

  class_<QDetectorConstruction, QDetectorConstruction*,
    bases<G4VUserDetectorConstruction> >
    ("QDetectorConstruction", "my detector")
    ;

  // ---
  def("Construct",  Construct);
}

