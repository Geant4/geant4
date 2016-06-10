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
// $Id: pyG4ApplicationState.cc 66892 2013-01-17 10:57:59Z gunter $
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
  enum_<G4ApplicationState>("G4ApplicationState")
    .value("G4State_PreInit",       G4State_PreInit)
    .value("G4State_Init",          G4State_Init)
    .value("G4State_Idle",          G4State_Idle)
    .value("G4State_GeomClosed",    G4State_GeomClosed)
    .value("G4State_EventProc",     G4State_EventProc)
    .value("G4State_Quit",          G4State_Quit)
    .value("G4State_Abort",         G4State_Abort)
    ;
}
