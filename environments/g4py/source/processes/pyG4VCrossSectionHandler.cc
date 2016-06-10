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
// $Id: pyG4VCrossSectionHandler.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4VCrossSectionHandler.cc
//
//                                         2008 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VCrossSectionHandler.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4Material.hh"
#include "G4Element.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4VCrossSectionHandler {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_Initialise, Initialise, 0, 8)

// FindValue
G4double (G4VCrossSectionHandler::*f1_FindValue)(G4int, G4double) const
  = &G4VCrossSectionHandler::FindValue;

G4double (G4VCrossSectionHandler::*f2_FindValue)(G4int, G4double, G4int) const
  = &G4VCrossSectionHandler::FindValue;

}

using namespace pyG4VCrossSectionHandler;

// ====================================================================
// module definition
// ====================================================================
void export_G4VCrossSectionHandler()
{
  class_<G4VCrossSectionHandler, boost::noncopyable>
    ("G4VCrossSectionHandler", "cross section handler", no_init)
    // ---
    .def("Initialise",           &G4VCrossSectionHandler::Initialise,
         f_Initialise())
    .def("SelectRandomElement",  &G4VCrossSectionHandler::SelectRandomElement,
         return_value_policy<reference_existing_object>())
    .def("SelectRandomShell",    &G4VCrossSectionHandler::SelectRandomShell)
    .def("FindValue",            f1_FindValue)
    .def("FindValue",            f2_FindValue)
    .def("ValueForMaterial",     &G4VCrossSectionHandler::ValueForMaterial)
    .def("LoadData",             &G4VCrossSectionHandler::LoadData)
    .def("LoadShellData",        &G4VCrossSectionHandler::LoadShellData)
    .def("PrintData",            &G4VCrossSectionHandler::PrintData)
    .def("Clear",                &G4VCrossSectionHandler::Clear)
    ;
}

