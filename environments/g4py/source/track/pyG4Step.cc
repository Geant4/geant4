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
// $Id: pyG4Step.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Step.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Step.hh"

using namespace boost::python;


// ====================================================================
// module definition
// ====================================================================
void export_G4Step()
{
  class_<G4Step, G4Step*>("G4Step", "step class")
    // ---
    .def("GetTrack",                 &G4Step::GetTrack,
         return_value_policy<reference_existing_object>())
    .def("GetPreStepPoint",          &G4Step::GetPreStepPoint,
         return_internal_reference<>())
    .def("GetPostStepPoint",         &G4Step::GetPostStepPoint,
         return_internal_reference<>())
    .def("GetTotalEnergyDeposit",    &G4Step::GetTotalEnergyDeposit)
    .def("GetStepLength",            &G4Step::GetStepLength)
    .def("GetDeltaPosition",         &G4Step::GetDeltaPosition)
    .def("GetDeltaTime",             &G4Step::GetDeltaTime)
    .def("GetDeltaMomentum",         &G4Step::GetDeltaMomentum)
    .def("GetDeltaEnergy",           &G4Step::GetDeltaEnergy)
    ;
}

