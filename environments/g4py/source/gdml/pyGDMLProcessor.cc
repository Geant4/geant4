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
// $Id: pyGDMLProcessor.cc,v 1.1 2007/11/13 10:03:23 kmura Exp $
// $Name: geant4-09-01 $
// ====================================================================
//   pyGDMLProcessor.cc
//
//   based on 2.10.0
//
//                                         2007 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Processor/GDMLProcessor.h"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4GDMLProcessor {

// GetWorldVolumeOfParsedFile
const G4VPhysicalVolume* (GDMLProcessor::*f1_GetWorldVolumeOfParsedFile)
  (const std::string&)= &GDMLProcessor::GetWorldVolumeOfParsedFile;

const G4VPhysicalVolume* (GDMLProcessor::*f2_GetWorldVolumeOfParsedFile)
  (const char*)= &GDMLProcessor::GetWorldVolumeOfParsedFile;

// GetWorldVolume
const G4VPhysicalVolume* (GDMLProcessor::*f1_GetWorldVolume)()
  = &GDMLProcessor::GetWorldVolume;

G4VPhysicalVolume* (GDMLProcessor::*f2_GetWorldVolume)(int)
  = &GDMLProcessor::GetWorldVolume;

}

using namespace pyG4GDMLProcessor;

// ====================================================================
// module definition
// ====================================================================
void export_GDMLProcessor()
{
  class_<GDMLProcessor, boost::noncopyable>
    ("GDMLProcessor", "GDML processor", no_init)
    // ---
    .def("GetInstance",  &GDMLProcessor::GetInstance,
         return_value_policy<reference_existing_object>())
    .staticmethod("GetInstance")
    // ---
    .def("GetWorldVolumeOfParsedFile",  f1_GetWorldVolumeOfParsedFile,
         return_value_policy<reference_existing_object>())
    .def("GetWorldVolumeOfParsedFile",  f2_GetWorldVolumeOfParsedFile,
         return_value_policy<reference_existing_object>())
    .def("GetWorldVolume",  f1_GetWorldVolume,
         return_value_policy<reference_existing_object>())
    .def("GetWorldVolume",  f2_GetWorldVolume,
         return_value_policy<reference_existing_object>())
    ;
}

