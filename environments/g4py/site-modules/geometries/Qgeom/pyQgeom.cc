//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: pyQgeom.cc,v 1.4 2006-06-29 15:29:29 gunter Exp $
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

