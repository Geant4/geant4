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
// $Id: pyG4ParticleDefinition.cc,v 1.6 2008-03-13 07:32:18 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ParticleDefinition.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
#include "G4ProcessManager.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4ParticleDefinition()
{
  class_<G4ParticleDefinition, G4ParticleDefinition*, boost::noncopyable> 
    ("G4ParticleDefinition", "particle definition", no_init)
    // ---
    .def("GetParticleName",    &G4ParticleDefinition::GetParticleName,
         return_value_policy<return_by_value>())
    .def("GetPDGMass",         &G4ParticleDefinition::GetPDGMass)
    .def("GetPDGWidth",        &G4ParticleDefinition::GetPDGWidth)
    .def("GetPDGCharge",       &G4ParticleDefinition::GetPDGCharge)
    .def("GetPDGSpin",         &G4ParticleDefinition::GetPDGSpin)
    .def("GetPDGiSpin",        &G4ParticleDefinition::GetPDGiSpin)
    .def("GetPDGiParity",      &G4ParticleDefinition::GetPDGiParity)
    .def("GetPDGiConjugation", &G4ParticleDefinition::GetPDGiConjugation)
    .def("GetPDGIsospin",      &G4ParticleDefinition::GetPDGIsospin)
    .def("GetPDGIsospin3",     &G4ParticleDefinition::GetPDGIsospin3)
    .def("GetPDGiIsospin",     &G4ParticleDefinition::GetPDGiIsospin)
    .def("GetPDGiIsospin3",    &G4ParticleDefinition::GetPDGiIsospin3)
    .def("GetPDGiGParity",     &G4ParticleDefinition::GetPDGiGParity)
    .def("GetParticleType",    &G4ParticleDefinition::GetParticleType,
         return_value_policy<return_by_value>())
    .def("GetParticleSubType", &G4ParticleDefinition::GetParticleSubType,
         return_value_policy<return_by_value>())
    .def("GetLeptonNumber",    &G4ParticleDefinition::GetLeptonNumber)
    .def("GetBaryonNumber",    &G4ParticleDefinition::GetBaryonNumber)
    .def("GetPDGEncoding",     &G4ParticleDefinition::GetPDGEncoding)
    .def("GetAntiPDGEncoding", &G4ParticleDefinition::GetAntiPDGEncoding)
    .def("GetQuarkContent",    &G4ParticleDefinition::GetQuarkContent) 
    .def("GetAntiQuarkContent",&G4ParticleDefinition::GetAntiQuarkContent)
    .def("IsShortLived",       &G4ParticleDefinition::IsShortLived)
    .def("GetPDGStable",       &G4ParticleDefinition::GetPDGStable)
    .def("SetPDGStable",       &G4ParticleDefinition::SetPDGStable)
    .def("GetPDGLifeTime",     &G4ParticleDefinition::GetPDGLifeTime)
    .def("SetPDGLifeTime",     &G4ParticleDefinition::SetPDGLifeTime)
    .def("GetDecayTable",      &G4ParticleDefinition::GetDecayTable,
         return_internal_reference<>())
    .def("SetDecayTable",      &G4ParticleDefinition::SetDecayTable)
    .def("GetProcessManager",  &G4ParticleDefinition::GetProcessManager,
         return_internal_reference<>())
    .def("SetProcessManager",  &G4ParticleDefinition::SetProcessManager)
    // cludge!! (G4ParticleTable object is sigleton!!)
    .def("GetParticleTable",   &G4ParticleDefinition::GetParticleTable,
         return_value_policy<reference_existing_object>()) 
    .def("DumpTable",          &G4ParticleDefinition::DumpTable)
#if G4VERSION_NUMBER >= 701
    .def("GetAtomicNumber",    &G4ParticleDefinition::GetAtomicNumber)
    .def("GetAtomicMass",      &G4ParticleDefinition::GetAtomicMass) 
#endif
    .def("SetVerboseLevel",    &G4ParticleDefinition::SetVerboseLevel)
    .def("GetVerboseLevel",    &G4ParticleDefinition::GetVerboseLevel)
    .def("SetApplyCutsFlag",   &G4ParticleDefinition::SetApplyCutsFlag)     
    .def("GetApplyCutsFlag",   &G4ParticleDefinition::GetApplyCutsFlag)     
    ;
}

