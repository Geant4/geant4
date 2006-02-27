// $Id: pyG4VProcess.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VProcess.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VProcess.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4VProcess {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_StorePhysicsTable, 
				       StorePhysicsTable, 2, 3);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_RetrievePhysicsTable, 
				       RetrievePhysicsTable, 2, 3);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetPhysicsTableFileName, 
				       GetPhysicsTableFileName, 3, 4);
};

using namespace pyG4VProcess;

// ====================================================================
// module definition
// ====================================================================
void export_G4VProcess()
{
  class_<G4VProcess, G4VProcess*, boost::noncopyable>
    ("G4VProcess", "base class for process", no_init)
    // ---
    // Note that only limited methods are exposed.
    .def("SetPILfactor",         &G4VProcess::SetPILfactor)
    .def("GetPILfactor",         &G4VProcess::GetPILfactor)
    .def("IsApplicable",         &G4VProcess::IsApplicable)
    .def("BuildPhysicsTable",    &G4VProcess::BuildPhysicsTable)
    .def("PreparePhysicsTable",  &G4VProcess::PreparePhysicsTable)
    .def("StorePhysicsTable",    &G4VProcess::StorePhysicsTable,
	 f_StorePhysicsTable())
    .def("RetrievePhysicsTable", &G4VProcess::RetrievePhysicsTable,
	 f_RetrievePhysicsTable())
    .def("GetPhysicsTableFileName", &G4VProcess::GetPhysicsTableFileName,
	 f_GetPhysicsTableFileName()
	 [return_value_policy<return_by_value>()])
    .def("GetProcessName",       &G4VProcess::GetProcessName,
	 return_value_policy<return_by_value>())  
    .def("GetProcessType",       &G4VProcess::GetProcessType)    
    .def("DumpInfo",             &G4VProcess::DumpInfo)
    .def("SetVerboseLevel",      &G4VProcess::SetVerboseLevel)
    .def("GetVerboseLevel",      &G4VProcess::GetVerboseLevel)
    ;
}
