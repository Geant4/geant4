// $Id: pyG4TouchableHistory.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VTouchable.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4TouchableHistory {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetVolume, GetVolume, 0, 1);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetTranslation,
				       GetTranslation, 0, 1);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetReplicaNumber, 
				       GetReplicaNumber, 0, 1);
}

using namespace pyG4TouchableHistory;

// ====================================================================
// module definition
// ====================================================================
void export_G4TouchableHistory()
{
  class_<G4TouchableHistory, G4TouchableHistory*, bases<G4VTouchable> >
    ("G4TouchableHistory", "touchable history class")
    // ---
    .def("GetVolume", &G4TouchableHistory::GetVolume,
	 f_GetVolume()
	 [return_value_policy<reference_existing_object>()])
    //.def(G4VSolid* GetSolid( G4int depth=0 ) const;
    .def("GetTranslation", &G4TouchableHistory::GetTranslation,
	 f_GetTranslation()
	 [return_value_policy<return_by_value>()])
    //.defG4RotationMatrix* GetRotation( G4int depth=0 ) const;
    .def("GetReplicaNumber", &G4TouchableHistory::GetReplicaNumber, 
	 f_GetReplicaNumber())
    .def("GetHistoryDepth",  &G4TouchableHistory::GetHistoryDepth)
    ;
}

