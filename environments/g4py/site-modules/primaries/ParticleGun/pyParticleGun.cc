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
// $Id: pyParticleGun.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyParticleGun.cc
//
//   [ParticleGun]
//   a site-module of Geant4Py
//
//   primary generator action with particle gun
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "ParticleGunAction.hh"
#include "G4ParticleGun.hh"
#include "G4RunManager.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyParticleGun {

G4ParticleGun* Construct()
{
  G4RunManager* runMgr= G4RunManager::GetRunManager();

  ParticleGunAction* pga= new ParticleGunAction;
  runMgr-> SetUserAction(pga);

  return pga-> GetParticleGun();
}

}

using namespace pyParticleGun;

// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(ParticleGun) {
  class_<ParticleGunAction, ParticleGunAction*,
    bases<G4VUserPrimaryGeneratorAction> >
    ("ParticleGunAction", "primary generator action with particle gun")
    .def("GetParticleGun", &ParticleGunAction::GetParticleGun,
         return_internal_reference<>())
    ;

  // ---
  def("Construct", Construct,
      return_value_policy<reference_existing_object>());

}

