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
// $Id: pyG4ProductionCutsTable.cc,v 1.4 2008-03-13 07:32:18 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ProductionCutsTable.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Material.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4ProductionCutsTable()
{
  class_<G4ProductionCutsTable, boost::noncopyable>
    ("G4ProductionCutsTable", "production cuts table", no_init)    
    .def("GetProductionCutsTable", 
         &G4ProductionCutsTable::GetProductionCutsTable,
         return_value_policy<reference_existing_object>())
    .staticmethod("GetProductionCutsTable")

    // internally used methods are limmitted to be exposed...

    // ---
    .def("GetLowEdgeEnergy",   &G4ProductionCutsTable::GetLowEdgeEnergy)
    .def("GetHighEdgeEnergy",  &G4ProductionCutsTable::GetHighEdgeEnergy)
    .def("SetEnergyRange",     &G4ProductionCutsTable::SetEnergyRange)
    .def("DumpCouples",        &G4ProductionCutsTable::DumpCouples)
    .def("IsModified",         &G4ProductionCutsTable::IsModified)
    // ---
#if G4VERSION_NUMBER >= 830
    .def("ConvertRangeToEnergy", &G4ProductionCutsTable::ConvertRangeToEnergy)
#endif
    // ---
    .def("SetVerboseLevel",    &G4ProductionCutsTable::SetVerboseLevel)
    .def("GetVerboseLevel",    &G4ProductionCutsTable::GetVerboseLevel)
    ;
}

