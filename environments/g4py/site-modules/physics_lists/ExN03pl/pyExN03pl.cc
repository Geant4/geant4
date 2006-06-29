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
// $Id: pyExN03pl.cc,v 1.3 2006-06-29 15:30:07 gunter Exp $
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
