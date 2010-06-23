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
// $Id: bindingEnergy.cc,v 1.9 2010-06-23 19:25:36 mkelsey Exp $
//
// 20100622 M. Kelsey -- Replace all functionally with call-through to
//		G4NucleiProperties.  Check for valid A/Z and return zero
//		without warning message.

#include "G4InuclSpecialFunctions.hh"
#include "G4NucleiProperties.hh"
#include "G4HadTmpUtil.hh"


G4double G4InuclSpecialFunctions::bindingEnergy(G4double A, G4double Z) {
  // NOTE:  Test condition copied from G4NucleiProperties.cc; not encapsulated
  if (A < 1 || Z < 0 || Z > A) return 0.;

  return G4NucleiProperties::GetBindingEnergy(G4lrint(A), G4lrint(Z));
}
