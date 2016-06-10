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
// $Id: pyExN03geom.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyExN03geom.cc
//
//   [ExN03geom]
//   a site-module of Geant4Py
//
//   geometry presented in ExN03 of Geant4 example
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4RunManager.hh"
#include "ExN03DetectorConstruction.hh"

using namespace boost::python;

typedef ExN03DetectorConstruction XXX;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyExN03geom {

void Construct()
{
  G4RunManager* runMgr= G4RunManager::GetRunManager();
  runMgr-> SetUserInitialization(new ExN03DetectorConstruction);
}

}

using namespace pyExN03geom;

// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(ExN03geom) {
  class_<ExN03DetectorConstruction, ExN03DetectorConstruction*,
    bases<G4VUserDetectorConstruction> >
    ("ExN03DetectorConstruction", "ExN03 detector")
    // ---
    .def("SetAbsorberMaterial",   &XXX::SetAbsorberMaterial)
    .def("SetAbsorberThickness",  &XXX::SetAbsorberThickness)
    .def("SetGapMaterial",        &XXX::SetGapMaterial)
    .def("SetGapThickness",       &XXX::SetGapThickness)
    .def("SetCalorSizeYZ",        &XXX::SetCalorSizeYZ)
    .def("SetNbOfLayers",         &XXX::SetNbOfLayers)
    .def("SetMagField",           &XXX::SetMagField)
    // ---
    .def("GetWorldSizeX",         &XXX::GetWorldSizeX)
    .def("GetWorldSizeYZ",        &XXX::GetWorldSizeYZ)
    .def("GetCalorThickness",     &XXX::GetCalorThickness)
    .def("GetCalorSizeYZ",        &XXX::GetCalorSizeYZ)
    .def("GetNbOfLayers",         &XXX::GetNbOfLayers)
    .def("GetAbsorberMaterial",   &XXX::GetAbsorberMaterial,
	 return_value_policy<reference_existing_object>())
    .def("GetAbsorberThickness",  &XXX::GetAbsorberThickness)
    .def("GetGapMaterial",        &XXX::GetGapMaterial,
	 return_value_policy<reference_existing_object>())
    .def("GetGapThickness",       &XXX::GetGapThickness)
    .def("GetphysiWorld",         &XXX::GetphysiWorld,
	 return_value_policy<reference_existing_object>())
    .def("GetAbsorber",           &XXX::GetAbsorber,
	 return_value_policy<reference_existing_object>())
    .def("GetGap",                &XXX::GetGap,
	 return_value_policy<reference_existing_object>())
    // ---
    .def("UpdateGeometry",        &XXX::UpdateGeometry)
    .def("PrintCalorParameters",  &XXX::PrintCalorParameters)
    ;

  // ---
  def("Construct",  Construct);
}

