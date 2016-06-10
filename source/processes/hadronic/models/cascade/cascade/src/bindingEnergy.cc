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
// $Id: bindingEnergy.cc 66241 2012-12-13 18:34:42Z gunter $
//
// 20100622 M. Kelsey -- Replace all functionally with call-through to
//		G4NucleiProperties.  Check for valid A/Z and return zero
//		without warning message.
// 20100914 M. Kelsey -- Migrate to integer A and Z

#include "G4InuclSpecialFunctions.hh"
#include "G4NucleiProperties.hh"


G4double G4InuclSpecialFunctions::bindingEnergy(G4int A, G4int Z) {
  // NOTE:  Test condition copied from G4NucleiProperties.cc; not encapsulated
  if (A < 1 || Z < 0 || Z > A) return 0.;

  return G4NucleiProperties::GetBindingEnergy(A, Z);
}
