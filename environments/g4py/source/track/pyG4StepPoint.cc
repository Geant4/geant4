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
// $Id: pyG4StepPoint.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
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
