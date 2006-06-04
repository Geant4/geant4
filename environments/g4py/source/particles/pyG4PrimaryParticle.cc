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
// $Id: pyG4PrimaryParticle.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4PrimaryParticle.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4PrimaryParticle.hh"
#include "G4ParticleDefinition.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4PrimaryParticle()
{
  class_<G4PrimaryParticle, G4PrimaryParticle*>
    ("G4PrimaryParticle", "primary particle")
    // ---
    .add_property("Px", &G4PrimaryParticle::GetPx)
    .add_property("Py", &G4PrimaryParticle::GetPy)
    .add_property("Pz", &G4PrimaryParticle::GetPz)
    // ---
    .def("Print",       &G4PrimaryParticle::Print)
    .def("GetPDGcode",  &G4PrimaryParticle::GetPDGcode)
    .def("GetG4code",   &G4PrimaryParticle::GetG4code,
         return_internal_reference<>())	 
    .def("GetMomentun", &G4PrimaryParticle::GetMomentum,
         return_value_policy<return_by_value>())
    .def("GetPx",       &G4PrimaryParticle::GetPx)
    .def("GetPy",       &G4PrimaryParticle::GetPy)
    .def("GetPz",       &G4PrimaryParticle::GetPz)
    .def("GetNext",     &G4PrimaryParticle::GetNext,
         return_internal_reference<>())
    .def("GetDaughter", &G4PrimaryParticle::GetNext,
         return_internal_reference<>())
    .def("GetTrackID",  &G4PrimaryParticle::GetTrackID)
    .def("GetMass",     &G4PrimaryParticle::GetMass)
    .def("GetCharge",   &G4PrimaryParticle::GetCharge)
    .def("GetPolarization", &G4PrimaryParticle::GetPolarization,
         return_value_policy<return_by_value>())
    .def("GetPolX",     &G4PrimaryParticle::GetPolX)
    .def("GetPolY",     &G4PrimaryParticle::GetPolY)
    .def("GetPolZ",     &G4PrimaryParticle::GetPolZ)
    .def("GetWeight",   &G4PrimaryParticle::GetWeight)
    .def("SetWeight",   &G4PrimaryParticle::SetWeight)
    .def("GetProperTime", &G4PrimaryParticle::GetProperTime)
    ;
}

