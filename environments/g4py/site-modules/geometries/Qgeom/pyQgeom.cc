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
// $Id: pyQgeom.cc,v 1.3 2006-06-04 21:36:35 kmura Exp $
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

