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
// Geant4 Version information

// Author: K.Murakami, 26.09.2005 - Created
// --------------------------------------------------------------------
#ifndef G4VERSION_HH
#define G4VERSION_HH 1

// Numbering rule for "G4VERSION_NUMBER":
// - The number is consecutive (i.e. 711) as an integer.
// - The meaning of each digit is as follows;
//
//   711
//   |--> major version number
//    |--> minor version number
//     |--> patch number

#ifndef G4VERSION_NUMBER
  #define G4VERSION_NUMBER 1111
#endif

#ifndef G4VERSION_TAG
  #define G4VERSION_TAG "$Name: geant4-11-01-patch-01 $"
#endif

// as variables

#include "G4String.hh"
#include "G4Types.hh"

#ifdef G4MULTITHREADED
static const G4String G4Version = "$Name: geant4-11-01-patch-01 [MT]$";
#else
static const G4String G4Version = "$Name: geant4-11-01-patch-01 $";
#endif
static const G4String G4Date = "(10-February-2023)";

#endif
