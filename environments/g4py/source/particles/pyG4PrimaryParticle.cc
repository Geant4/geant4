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
// $Id: pyG4PrimaryParticle.cc,v 1.4 2006-06-29 15:34:40 gunter Exp $
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

