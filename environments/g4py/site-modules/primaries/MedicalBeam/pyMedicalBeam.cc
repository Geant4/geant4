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
// $Id: pyMedicalBeam.cc,v 1.4 2006-06-29 15:30:14 gunter Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyMedicalBeam.cc
//
//   [MedicalBeam]
//   a site-module of Geant4Py
//
//   primary generator action for medical beam
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "MedicalBeam.hh"
#include "G4ParticleTable.hh"
#include "G4RunManager.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyMedicalBeam {

MedicalBeam* Construct()
{
  G4RunManager* runMgr= G4RunManager::GetRunManager();

  MedicalBeam* medicalbeam= new MedicalBeam();
  runMgr-> SetUserAction(medicalbeam);

  return medicalbeam;
}


///////////////////////////////////////////////////////////////////
void SetParticleByName(MedicalBeam* beam, const std::string& pname)
///////////////////////////////////////////////////////////////////
{
  G4ParticleTable* particleTable= G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* pd= particleTable-> FindParticle(pname);
  if (pd != 0) {
    beam-> SetParticleDefinition(pd);
  } else {
    G4cout << "*** \"" << pname << "\" is not registered "
           << "in available particle list" << G4endl;
  }
}

////////////////////////////////////////////////
std::string GetParticleByName(MedicalBeam* beam)
////////////////////////////////////////////////
{
  const G4ParticleDefinition* pd= beam-> GetParticleDefinition();

  if(pd==0) return std::string("None");
  else return (pd-> GetParticleName()).c_str();
}

////////////////////////////////////////////////////////
void f_SetFieldXY(MedicalBeam* beam, const list& listXY)
////////////////////////////////////////////////////////
{
  G4double fx= extract<double>(listXY[0]);
  G4double fy= extract<double>(listXY[1]);
  beam-> SetFieldXY(fx, fy);
}


////////////////////////////////////
list f_GetFieldXY(MedicalBeam* beam)
////////////////////////////////////
{
  list listFieldXY;

  listFieldXY.append(beam-> GetFieldX());
  listFieldXY.append(beam-> GetFieldY());

  return listFieldXY;
}

};

using namespace pyMedicalBeam;

// ====================================================================
//   Expose to Python
// ====================================================================

BOOST_PYTHON_MODULE(MedicalBeam) {
  class_<MedicalBeam, MedicalBeam*,
    bases<G4VUserPrimaryGeneratorAction> >
    ("MedicalBeam", "primary generator action with medical beam")
    // ---
    .add_property("particle", GetParticleByName, SetParticleByName)
    .def("SetParticleByName", SetParticleByName)
    .def("GetParticleByName", GetParticleByName)
    // ---
    .add_property("kineticE", &MedicalBeam::GetKineticE, 
		              &MedicalBeam::SetKineticE)
    .def("SetKineticE",       &MedicalBeam::SetKineticE)
    .def("GetKineticE",       &MedicalBeam::GetKineticE)
    // ---
    .add_property("sourcePosition", &MedicalBeam::GetSourcePosition,
                  		    &MedicalBeam::SetSourcePosition)
    .def("SetSourcePosition",       &MedicalBeam::SetSourcePosition)
    .def("GetSourcePosition",       &MedicalBeam::GetSourcePosition)
    // ---
    .add_property("fieldShape", &MedicalBeam::GetFieldShape,
		                &MedicalBeam::SetFieldShape)
    .def("SetFieldShape",       &MedicalBeam::SetFieldShape)
    .def("GetFieldShape",       &MedicalBeam::GetFieldShape)
    // ---
    .add_property("SSD", &MedicalBeam::GetSSD, &MedicalBeam::SetSSD)
    .def("SetSSD",       &MedicalBeam::SetSSD)
    .def("GetSSD",       &MedicalBeam::GetSSD)
    // ----
    .add_property("fieldXY", f_GetFieldXY, f_SetFieldXY)
    .def("SetFieldXY",       f_SetFieldXY)
    .def("GetFieldXY",       f_GetFieldXY)
    .def("GetFieldX",        &MedicalBeam::GetFieldX)
    .def("GetFieldY",        &MedicalBeam::GetFieldY)
    // ---
    .add_property("fieldR", &MedicalBeam::GetFieldR, &MedicalBeam::SetFieldR)
    .def("SetFieldR",       &MedicalBeam::SetFieldR)
    .def("GetFieldR",       &MedicalBeam::GetFieldR)
    ;
 
  // enums...
  enum_<MedicalBeam::FieldShape>("FieldShape")
    .value("SQUARE", MedicalBeam::SQUARE)
    .value("CIRCLE", MedicalBeam::CIRCLE)
    ;
  
  // ---
  def("Construct",  Construct,
      return_value_policy<reference_existing_object>());
}

