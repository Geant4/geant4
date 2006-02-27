// $Id: pyG4StepStatus.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4StepStatus.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4StepStatus.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4StepStatus()
{
  enum_<G4StepStatus>("G4StepStatus")
    .value("fWorldBoundary",         fWorldBoundary)
    .value("fGeomBoundary",          fGeomBoundary)
    .value("fAtRestDoItProc",        fAtRestDoItProc)
    .value("fAlongStepDoItProc",     fAlongStepDoItProc)
    .value("fPostStepDoItProc",      fPostStepDoItProc)
    .value("fUserDefinedLimit",      fUserDefinedLimit)
    .value("fExclusivelyForcedProc", fExclusivelyForcedProc)
    .value("fUndefined",             fUndefined)
    ;
}

