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
// $Id: pyG4UIparameter.cc 66892 2013-01-17 10:57:59Z gunter $
// ====================================================================
//   pyG4UIparameter.cc
//
//                                         2006 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4UIparameter.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4UIparameter()
{
  class_<G4UIparameter, G4UIparameter*>
    ("G4UIparameter", "UI parameter")
    // constructors
    .def(init<char>())
    .def(init<const char*, char, G4bool>())
    // ---
    .def("List",                   &G4UIparameter::List)
    .def("GetDefaultValue",        &G4UIparameter::GetDefaultValue)
    .def("GetParameterType",       &G4UIparameter::GetParameterType)
    .def("GetParameterRange",      &G4UIparameter::GetParameterRange)
    .def("GetParameterName",       &G4UIparameter::GetParameterName)
    .def("GetParameterCandidates", &G4UIparameter::GetParameterCandidates)
    .def("IsOmittable",            &G4UIparameter::IsOmittable)
    .def("GetCurrentAsDefault",    &G4UIparameter::GetCurrentAsDefault)
    .def("GetParameterGuidance",   &G4UIparameter::GetParameterGuidance)
    ;
}
