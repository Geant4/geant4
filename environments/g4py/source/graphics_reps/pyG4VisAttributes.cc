// $Id: pyG4VisAttributes.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VisAttributes.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VisAttributes.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4VisAttributes {

// SetColor()
void(G4VisAttributes::*f1_SetColor)(const G4Color&)= 
  &G4VisAttributes::SetColor;

void(G4VisAttributes::*f2_SetColor)(G4double, G4double, G4double, G4double)= 
  &G4VisAttributes::SetColor;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_SetColor, SetColor, 3, 4);

// SetColour()
void(G4VisAttributes::*f1_SetColour)(const G4Colour&)= 
  &G4VisAttributes::SetColour;

void(G4VisAttributes::*f2_SetColour)(G4double, G4double, G4double, G4double)= 
  &G4VisAttributes::SetColour;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_SetColour, SetColour, 3, 4);

};

using namespace pyG4VisAttributes;

// ====================================================================
// module definition
// ====================================================================
void export_G4VisAttributes()
{
  scope in_G4VisAttributes= 
    class_<G4VisAttributes, G4VisAttributes*>
    ("G4VisAttributes", "visualization attributes")
    // constructors
    .def(init<G4bool>())
    .def(init<const G4Colour&>())
    .def(init<G4bool, const G4Colour&>())
    // ---
    .def_readonly("Invisible",     &G4VisAttributes::Invisible)
    .def("GetInvisible",           &G4VisAttributes::GetInvisible,
	 return_value_policy<reference_existing_object>())
    .staticmethod("GetInvisible")    
    .def("IsVisible",              &G4VisAttributes::IsVisible)
    .def("IsDaughtersInvisible",   &G4VisAttributes::IsDaughtersInvisible)
    // ---
    .def("GetColour",              &G4VisAttributes::GetColour,
         return_internal_reference<>())
    .def("GetColor",               &G4VisAttributes::GetColor,
         return_internal_reference<>())
    // ---
    .def("GetLineStyle",           &G4VisAttributes::GetLineStyle)
    .def("GetLineWidth",           &G4VisAttributes::GetLineWidth)
    .def("IsForceDrawingStyle",    &G4VisAttributes::IsForceDrawingStyle)
    .def("GetForcedDrawingStyle",  &G4VisAttributes::GetForcedDrawingStyle)
    .def("IsForceAuxEdgeVisible",  &G4VisAttributes::IsForceAuxEdgeVisible)
    .def("SetVisibility",          &G4VisAttributes::SetVisibility)
    .def("SetDaughtersInvisible",  &G4VisAttributes::SetDaughtersInvisible)
    .def("SetColor",               f1_SetColor)
    .def("SetColour",              f1_SetColour)
    .def("SetColor",               f2_SetColor,      f_SetColor())
    .def("SetColour",              f2_SetColour,     f_SetColour())
    .def("SetLineStyle",           &G4VisAttributes::SetLineStyle)
    .def("SetLineWidth",           &G4VisAttributes::SetLineWidth)
    .def("SetForceWireframe",      &G4VisAttributes::SetForceWireframe)
    .def("SetForceSolid",          &G4VisAttributes::SetForceSolid)
    .def("SetForceAuxEdgeVisible", &G4VisAttributes::SetForceAuxEdgeVisible)
    .def("SetAttValues",           &G4VisAttributes::SetAttValues)
    .def("SetAttDefs",             &G4VisAttributes::SetAttDefs)
    // operators
    .def(self_ns::str(self))
    .def(self == self)
    .def(self != self)
    ;

  // enum LineStyle
  enum_<G4VisAttributes::LineStyle>("LineStyle")
    .value("unbroken",   G4VisAttributes::unbroken)
    .value("dashed",     G4VisAttributes::dashed)
    .value("dotted",     G4VisAttributes::dotted)
    ;

  // enum ForcedDrawingStyle
  enum_<G4VisAttributes::ForcedDrawingStyle>("ForcedDrawingStyle")
    .value("wireframe",  G4VisAttributes::wireframe)
    .value("solid",      G4VisAttributes::solid)
    ;

}
