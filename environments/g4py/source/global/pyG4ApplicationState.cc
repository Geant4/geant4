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
// $Id: pyG4ApplicationState.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
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
