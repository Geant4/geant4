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
// $Id: pyExN03pl.cc,v 1.2 2006-06-04 21:36:35 kmura Exp $
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
