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
// $Id: pyG4ParticleGun.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4ParticleGun.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Event.hh"

using namespace boost::python;

#if G4VERSION_NUMBER < 910
// ====================================================================
// miscs
// ====================================================================
// What a hell!

////////////////////////////////////////////////////////
G4ParticleGun::G4ParticleGun(const G4ParticleGun &right)
////////////////////////////////////////////////////////
{
  *this= right;
}

/////////////////////////////////////////////////////////////////////////
const G4ParticleGun& G4ParticleGun::operator=(const G4ParticleGun &right)
/////////////////////////////////////////////////////////////////////////
{
   NumberOfParticlesToBeGenerated= right.NumberOfParticlesToBeGenerated;
   particle_definition= right.particle_definition;
   particle_momentum_direction= right.particle_momentum_direction;
   particle_energy= right.particle_energy;
   particle_charge= right.particle_charge;
   particle_polarization= right.particle_polarization;

   return *this;
}

/////////////////////////////////////////////////////////////////
G4int G4ParticleGun::operator==(const G4ParticleGun &right) const
/////////////////////////////////////////////////////////////////
{
  return 0;
}

/////////////////////////////////////////////////////////////////
G4int G4ParticleGun::operator!=(const G4ParticleGun &right) const
/////////////////////////////////////////////////////////////////
{
  return 0;
}

#endif


// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4ParticleGun {

#if G4VERSION_NUMBER >= 910
// SetParticleMomentum
void (G4ParticleGun::*f1_SetParticleMomentum)(G4double)
  = &G4ParticleGun::SetParticleMomentum;
void (G4ParticleGun::*f2_SetParticleMomentum)(G4ParticleMomentum)
  = &G4ParticleGun::SetParticleMomentum;
#endif


////////////////////////////////////////////////////////////////////
void SetParticleByName(G4ParticleGun* gun, const std::string& pname)
////////////////////////////////////////////////////////////////////
{
  G4ParticleTable* particleTable= G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* pd= particleTable-> FindParticle(pname);
  if (pd != 0) {
    gun-> SetParticleDefinition(pd);
  } else {
    G4cout << "*** \"" << pname << "\" is not registered "
	   << "in available particle list" << G4endl;
  }
}

/////////////////////////////////////////////////
std::string GetParticleByName(G4ParticleGun* gun)
/////////////////////////////////////////////////
{
  const G4ParticleDefinition* pd= gun-> GetParticleDefinition();
  return (pd-> GetParticleName()).c_str();
}

}

using namespace pyG4ParticleGun;

// ====================================================================
// module definition
// ====================================================================
void export_G4ParticleGun()
{
#if G4VERSION_NUMBER < 910
  class_<G4ParticleGun>
#else
    class_<G4ParticleGun, boost::noncopyable>
#endif
    ("G4ParticleGun", "particle gun")
    // constructor
    .def(init<G4int>())
    .def(init<G4ParticleDefinition*>())
    .def(init<G4ParticleDefinition*, G4int>())
    // ---
    .def("GeneratePrimaryVertex", &G4ParticleGun::GeneratePrimaryVertex)
    .def("SetParticleDefinition", &G4ParticleGun::SetParticleDefinition)
    .def("GetParticleDefinition", &G4ParticleGun::GetParticleDefinition,
    	 return_value_policy<reference_existing_object>())
#if G4VERSION_NUMBER >= 910
    .def("SetParticleMomentum",   f1_SetParticleMomentum)
    .def("SetParticleMomentum",   f2_SetParticleMomentum)
#else
    .def("SetParticleMomentum",   &G4ParticleGun::SetParticleMomentum)
#endif
    .def("SetParticleMomentumDirection",
	 &G4ParticleGun::SetParticleMomentumDirection)
    .def("GetParticleMomentumDirection",
	 &G4ParticleGun::GetParticleMomentumDirection)
    .def("SetParticleEnergy",     &G4ParticleGun::SetParticleEnergy)
    .def("GetParticleEnergy",     &G4ParticleGun::GetParticleEnergy)
    .def("SetParticleCharge",     &G4ParticleGun::SetParticleCharge)
    .def("GetParticleCharge",     &G4ParticleGun::GetParticleCharge)
    .def("SetParticlePolarization", &G4ParticleGun::SetParticlePolarization)
    .def("GetParticlePolarization", &G4ParticleGun::GetParticlePolarization)
    .def("SetNumberOfParticles",  &G4ParticleGun::SetNumberOfParticles)
    .def("GetNumberOfParticles",  &G4ParticleGun::GetNumberOfParticles)
    .def("SetParticlePosition",   &G4ParticleGun::SetParticlePosition)
    .def("GetParticlePosition",   &G4ParticleGun::GetParticlePosition)
    .def("SetParticleTime",       &G4ParticleGun::SetParticleTime)
    .def("GetParticleTime",       &G4ParticleGun::GetParticleTime)
    .def("SetParticleByName",     SetParticleByName)
    .def("GetParticleByName",     GetParticleByName)
    ;
}
