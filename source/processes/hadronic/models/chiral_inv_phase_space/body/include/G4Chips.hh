// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Chips.hh,v 1.2 1999-12-15 14:52:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4Chips_h
#define G4Chips_h 1

// ------------------------------------------------------------
//      GEANT 4 namespace header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4Chips ----------------
//          by Mikhail Kossov, September 1999.
//      namespace for constants of the CHIPS Model
// ------------------------------------------------------------

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
