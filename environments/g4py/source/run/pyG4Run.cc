// $Id: pyG4Run.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
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
