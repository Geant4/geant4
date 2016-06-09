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
// $Id: pyG4DynamicParticle.cc,v 1.6 2008-03-13 07:32:18 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4DynamicParticle.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#include "G4DynamicParticle.hh"
#if G4VERSION_NUMBER <= 701
#include "G4PrimaryParticle.hh"
#endif

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4DynamicParticle()
{
  class_<G4DynamicParticle, G4DynamicParticle*>
    ("G4DynamicParticle", "dynamic particle")
    // ---
    .def("GetMomentumDirection", &G4DynamicParticle::GetMomentumDirection,
	 return_value_policy<return_by_value>())
    .def("GetMomentum",          &G4DynamicParticle::GetMomentum,
	 return_value_policy<return_by_value>())
    //.def("Get4Momentum",       &G4DynamicParticle::Get4Momentum,
    //return_value_policy<return_by_value>())
    .def("GetTotalMomentum",     &G4DynamicParticle::GetTotalMomentum)
    .def("GetTotalEnergy",       &G4DynamicParticle::GetTotalEnergy)
    .def("GetKineticEnergy",     &G4DynamicParticle::GetKineticEnergy)
    .def("GetProperTime",        &G4DynamicParticle::GetProperTime)
    .def("GetPolarization",      &G4DynamicParticle::GetPolarization,
	 return_value_policy<return_by_value>())
    .def("GetMass",              &G4DynamicParticle::GetMass)
    .def("GetCharge",            &G4DynamicParticle::GetCharge)
    //.def("GetElectronOccupancy", &G4DynamicParticle::GetElectronOccupancy,
    //return_internal_reference<>())
    .def("GetTotalOccupancy",    &G4DynamicParticle::GetTotalOccupancy)
    .def("GetOccupancy",         &G4DynamicParticle::GetOccupancy)
    .def("GetDefinition",        &G4DynamicParticle::GetDefinition,
	 return_internal_reference<>())
    .def("GetPreAssignedDecayProperTime", 
	 &G4DynamicParticle::GetPreAssignedDecayProperTime)
    .def("DumpInfo",             &G4DynamicParticle::DumpInfo)
    .def("SetVerboseLevel",      &G4DynamicParticle::SetVerboseLevel)
    .def("GetVerboseLevel",      &G4DynamicParticle::GetVerboseLevel)
    .def("GetPrimaryParticle",   &G4DynamicParticle::GetPrimaryParticle,
	 return_internal_reference<>())
#if G4VERSION_NUMBER >= 710
    .def("GetPDGcode",           &G4DynamicParticle::GetPDGcode)
#endif
    ;   
}

