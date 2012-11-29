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
//
// $Id$
//
//      ---------------- G4Chips ----------------
//          by Mikhail Kossov, September 1999.
//      namespace for constants of the CHIPS Model
// ------------------------------------------------------------

#ifndef G4Chips_h
#define G4Chips_h 1

// >->-> D O E S   N O T   W O R K   -   N O T   U S E D    N O W <-<-<
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
