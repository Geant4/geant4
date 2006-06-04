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
// $Id: pyG4Run.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Run.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Run.hh"
#include "G4HCtable.hh"
#include "G4DCtable.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Run()
{
  class_<G4Run, G4Run*>("G4Run", "run class")
    // ---
    .def("GetRunID",         &G4Run::GetRunID)
    .def("SetRunID",         &G4Run::SetRunID)
    .def("GetNumberOfEvent", &G4Run::GetNumberOfEvent)
    .def("GetNumberOfEventToBeProcessed", 
	 &G4Run::GetNumberOfEventToBeProcessed)
    .def("SetNumberOfEventToBeProcessed", 
	 &G4Run::SetNumberOfEventToBeProcessed)
    ;

    // reduced functionality...
    //.def("RecordEvent",      &G4Run::RecordEvent) // virtual
    //.def("GetHCtable",       &G4Run::GetHCtable,
    //return_internal_reference<>())
    //.def("SetHCtable",       &G4Run::SetHCtable)
    //.def("GetDCtable",       &G4Run::GetDCtable,
    //return_internal_reference<>())
    //.def("SetDCtable",       &G4Run::SetDCtable)	 

}
