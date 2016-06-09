//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4PyCoutDestination.cc,v 1.3 2006-06-29 15:33:03 gunter Exp $
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


////////////////////////////////////////////////////////////////////
G4int G4PyCoutDestination::ReceiveG4cout(const G4String& coutString)
////////////////////////////////////////////////////////////////////
{
  PySys_WriteStdout("%s", coutString.c_str());
  return 0;
}


////////////////////////////////////////////////////////////////////
G4int G4PyCoutDestination::ReceiveG4cerr(const G4String& cerrString)
////////////////////////////////////////////////////////////////////
{
  PySys_WriteStderr("%s", cerrString.c_str());
  return 0;
}

