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
// $Id: pyParticleGun.cc,v 1.3 2006-06-04 21:36:35 kmura Exp $
// $Name: not supported by cvs2svn $
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

};

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

