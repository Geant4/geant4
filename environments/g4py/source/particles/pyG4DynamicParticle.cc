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
// $Id: pyG4DynamicParticle.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4DynamicParticle.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4DynamicParticle.hh"
#include "pyG4Version.hh"

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

