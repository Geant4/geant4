// $Id: pyG4ApplicationState.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ApplicationState.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4ApplicationState.hh"

using namespace boost::python;

// ====================================================================
//
// class description
//
// ====================================================================
void export_G4ApplicationState()
{
  enum_<G4ApplicationState>("G4APplicationState")
    .value("G4State_PreInit",       G4State_PreInit)
    .value("G4State_Init",          G4State_Init)
    .value("G4State_Idle",          G4State_Idle)
    .value("G4State_GeomClosed",    G4State_GeomClosed)
    .value("G4State_EventProc",     G4State_EventProc)
    .value("G4State_Quit",          G4State_Quit)
    .value("G4State_Abort",         G4State_Abort)
    ;
}
