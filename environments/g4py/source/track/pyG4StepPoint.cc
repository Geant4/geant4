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
// $Id: pyG4StepPoint.cc 66892 2013-01-17 10:57:59Z gunter $
// ====================================================================
//   pyG4StepPoint.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4StepPoint.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4StepPoint()
{
  class_<G4StepPoint, G4StepPoint*>("G4StepPoint", "step point class")
    // ---
    .def("GetPosition",           &G4StepPoint::GetPosition,
	 return_value_policy<return_by_value>())
    .def("GetLocalTime",          &G4StepPoint::GetLocalTime)
    .def("GetGlobalTime",         &G4StepPoint::GetGlobalTime)
    .def("GetProperTime",         &G4StepPoint::GetProperTime)
    .def("GetMomentumDirection",  &G4StepPoint::GetMomentumDirection,
	 return_value_policy<return_by_value>())
    .def("GetMomentum",           &G4StepPoint::GetMomentum,
	 return_value_policy<return_by_value>())
    .def("GetTotalEnergy",        &G4StepPoint::GetTotalEnergy)
    .def("GetKineticEnergy",      &G4StepPoint::GetKineticEnergy)
    .def("GetVelocity",           &G4StepPoint::GetVelocity)
    .def("GetBeta",               &G4StepPoint::GetBeta)
    .def("GetGamma",              &G4StepPoint::GetGamma)
    //.def("GetTouchable",          &G4StepPoint::GetTouchable)
    //.def("GetMaterial",           &G4StepPoint::GetMaterial)
    .def("GetPolarization",       &G4StepPoint::GetPolarization,
	 return_value_policy<return_by_value>())
    .def("GetStepStatus",         &G4StepPoint::GetStepStatus)
    //.def("GetProcessDefinedStep", &G4StepPoint::GetProcessDefinedStep)
    .def("GetMass",               &G4StepPoint::GetMass)
    .def("GetCharge",             &G4StepPoint::GetCharge)
    .def("GetWeight",             &G4StepPoint::GetWeight)
    ;
}
