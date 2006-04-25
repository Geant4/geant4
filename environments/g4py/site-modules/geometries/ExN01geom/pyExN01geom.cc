// $Id: pyExN01geom.cc,v 1.2 2006-04-25 10:25:12 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyExN01geom.cc
//
//   [Qgeom]
//   a site-module of Geant4Py
//
//   geometry presented in ExN01 of Geant4 example
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4RunManager.hh"
#include "ExN01DetectorConstruction.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyExN03geom {

void Construct()
{
  G4RunManager* runMgr= G4RunManager::GetRunManager();
  runMgr-> SetUserInitialization(new ExN01DetectorConstruction);
}

};

using namespace pyExN03geom;

// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(ExN01geom) {

  class_<ExN01DetectorConstruction, ExN01DetectorConstruction*,
    bases<G4VUserDetectorConstruction> >
    ("ExN01DetectorConstruction", "ExN01 detector")
    ;

  // ---
  def("Construct",  Construct);
}

