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
// $Id: pyG4LossTableManager.cc,v 1.5 2008-03-13 07:32:18 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4LossTableManager.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#include "G4LossTableManager.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4LossTableManager()
{
  class_<G4LossTableManager, boost::noncopyable>
    ("G4LossTableManager", "energy loss table manager", no_init)
    .def("Instance",    &G4LossTableManager::Instance,
         return_value_policy<reference_existing_object>())
    .staticmethod("Instance")

    // internally used methods are limmitted to be exposed...

    // ----
    .def("SetLossFluctuations",  &G4LossTableManager::SetLossFluctuations)
    .def("SetSubCutoff",         &G4LossTableManager::SetSubCutoff)
    .def("SetIntegral",          &G4LossTableManager::SetIntegral)
    .def("SetRandomStep",        &G4LossTableManager::SetRandomStep)
    .def("SetMinSubRange",       &G4LossTableManager::SetMinSubRange)
    .def("SetMinEnergy",         &G4LossTableManager::SetMinEnergy)
    .def("SetMaxEnergy",         &G4LossTableManager::SetMaxEnergy)

#if G4VERSION_NUMBER >= 810
    .def("SetMaxEnergyForCSDARange", 
         &G4LossTableManager::SetMaxEnergyForCSDARange)
#else
    .def("SetMaxEnergyForPreciseRange", 
         &G4LossTableManager::SetMaxEnergyForPreciseRange)
#endif

    .def("SetMaxEnergyForMuons", &G4LossTableManager::SetMaxEnergyForMuons)

#if G4VERSION_NUMBER >= 810
    .def("SetStepFunction",      &G4LossTableManager::SetStepFunction)
    .def("SetBuildCSDARange",    &G4LossTableManager::SetBuildCSDARange)
    .def("SetLPMFlag",           &G4LossTableManager::SetLPMFlag)
    .def("SetBremsstrahlungTh",  &G4LossTableManager::SetBremsstrahlungTh)
#else
    .def("SetStepLimits",        &G4LossTableManager::SetStepLimits)
    .def("SetBuildPreciseRange", &G4LossTableManager::SetBuildPreciseRange)
#endif
    .def("SetVerbose",           &G4LossTableManager::SetVerbose)

    // ---
#if G4VERSION_NUMBER >= 810
    .def("BuildCSDARange",       &G4LossTableManager::BuildCSDARange)
    .def("LPMFlag",              &G4LossTableManager::LPMFlag)
    .def("BremsstrahlungTh",     &G4LossTableManager::BremsstrahlungTh)

#if G4VERSION_NUMBER <= 830
    .def("SetMscStepLimitation", &G4LossTableManager::SetMscStepLimitation)
    .def("MscFlag",              &G4LossTableManager::MscFlag)
    .def("FacRange",             &G4LossTableManager::FacRange)
#endif

#else
    .def("BuildPreciseRange",    &G4LossTableManager::BuildPreciseRange)
#endif

    // ---
#if G4VERSION_NUMBER == 830
    .def("SetMscLateralDisplacement", 
         &G4LossTableManager::SetMscLateralDisplacement)
    .def("SetSkin",              &G4LossTableManager::SetSkin)
#endif

#if G4VERSION_NUMBER >= 830
    .def("SetLinearLossLimit",   &G4LossTableManager::SetLinearLossLimit)
#endif

    ;

}

