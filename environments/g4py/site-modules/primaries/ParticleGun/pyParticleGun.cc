// $Id: pyParticleGun.cc,v 1.1 2006-02-27 09:52:54 kmura Exp $
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

