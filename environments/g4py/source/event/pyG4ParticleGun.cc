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
// $Id: pyG4ParticleGun.cc,v 1.3 2006-06-04 21:34:28 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ParticleGun.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Event.hh"

using namespace boost::python;

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

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4ParticleGun {

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

};

using namespace pyG4ParticleGun;

// ====================================================================
// module definition
// ====================================================================
void export_G4ParticleGun()
{
  class_<G4ParticleGun>
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
    .def("SetParticleMomentum",   &G4ParticleGun::SetParticleMomentum)
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
