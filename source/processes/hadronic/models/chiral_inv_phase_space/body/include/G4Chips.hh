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
//
// $Id: G4Chips.hh,v 1.9 2001-11-21 11:47:21 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4Chips ----------------
//          by Mikhail Kossov, September 1999.
//      namespace for constants of the CHIPS Model
// ------------------------------------------------------------

#ifndef G4Chips_h
#define G4Chips_h 1

// >>> D O E S   N O T   W O R K   -   N O T   U S E D    N O W <<<
// Instead the "static const G4double" are used in the member functions

namespace G4Chips
{
  // To have 3 quarks in Nucleon Temperature should be < M_N/4 (234 MeV)
  const G4double  Temperature = 200.;  // Temperature of Quasmon (constant of the model)
}

// Now the same constants or types should be made global (if it is necessary)
//
const double& G4Chips::Temperature Temperature;
// typedef G4Chips::NewType Newtype;
