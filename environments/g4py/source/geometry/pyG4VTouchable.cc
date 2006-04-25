// $Id: pyG4VTouchable.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
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

