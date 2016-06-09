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
// $Id: pyG4VTouchable.cc,v 1.4 2006-06-29 15:32:57 gunter Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VTouchable.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4VTouchable {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetTranslation,
				       GetTranslation, 0, 1);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetRotation,
				       GetRotation, 0, 1);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetVolume, GetVolume, 0, 1);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetSolid, GetSolid, 0, 1);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetReplicaNumber, 
				       GetReplicaNumber, 0, 1);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_MoveUpHistory, MoveUpHistory, 0, 1);
  
}

using namespace pyG4VTouchable;

// ====================================================================
// module definition
// ====================================================================
void export_G4VTouchable()
{
  class_<G4VTouchable, G4VTouchable*, boost::noncopyable>
    ("G4VTouchable", "touchable class", no_init)
    .def("GetTranslation", &G4VTouchable::GetTranslation,
	 f_GetTranslation()
	 [return_value_policy<return_by_value>()])
    .def("GetRotation", &G4VTouchable::GetRotation,
	 f_GetRotation()
	 [return_value_policy<reference_existing_object>()])
    .def("GetVolume", &G4VTouchable::GetVolume,
	 f_GetVolume()
	 [return_value_policy<reference_existing_object>()])
    .def("GetSolid", &G4VTouchable::GetSolid,
	 f_GetSolid()
	 [return_value_policy<reference_existing_object>()])
    .def("GetReplicaNumber", &G4VTouchable::GetReplicaNumber, 
	 f_GetReplicaNumber())
    .def("GetHistoryDepth",  &G4VTouchable::GetHistoryDepth)
    .def("MoveUpHistory", &G4VTouchable::MoveUpHistory, 
	 f_MoveUpHistory())
    ;
}

