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
// $Id: G4PyCoutDestination.cc,v 1.2 2006-06-04 21:34:29 kmura Exp $
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

