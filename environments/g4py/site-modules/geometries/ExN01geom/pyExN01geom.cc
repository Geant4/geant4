//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: pyExN01geom.cc,v 1.3 2006-06-04 21:36:34 kmura Exp $
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

