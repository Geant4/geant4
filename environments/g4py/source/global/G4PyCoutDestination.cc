// $Id: G4PyCoutDestination.cc,v 1.1 2006-04-25 08:13:51 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   G4PyCoutDistination.cc
//
//                                         2006 Q
// ====================================================================
#include <Python.h>
#include "G4PyCoutDestination.hh"

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////////////////////
G4PyCoutDestination::G4PyCoutDestination()
//////////////////////////////////////////
{
}

///////////////////////////////////////////
G4PyCoutDestination::~G4PyCoutDestination()
///////////////////////////////////////////
{
}


/////////////////////////////////////////////////////////////
G4int G4PyCoutDestination::ReceiveG4cout(G4String coutString)
/////////////////////////////////////////////////////////////
{
  PySys_WriteStdout("%s", coutString.c_str());
  return 0;
}


/////////////////////////////////////////////////////////////
G4int G4PyCoutDestination::ReceiveG4cerr(G4String cerrString)
/////////////////////////////////////////////////////////////
{
  PySys_WriteStderr("%s", cerrString.c_str());
  return 0;
}

