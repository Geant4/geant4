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
// 20110922  M. Kelsey -- Created to hold implementation of stream operator<<.
// 20121009  M. Kelsey -- Reduce output to single line, with no EOL.
// 20130622  M. Kelsey -- Add initialization from G4Fragment

#include "G4ExitonConfiguration.hh"
#include "G4Fragment.hh"
#include "G4ios.hh"


// Get information from G4Fragment

void G4ExitonConfiguration::fill(const G4Fragment& frag) {
  protonQuasiParticles  = frag.GetNumberOfCharged();
  neutronQuasiParticles = frag.GetNumberOfParticles() - protonQuasiParticles;
  protonHoles  = frag.GetNumberOfChargedHoles();
  neutronHoles = frag.GetNumberOfHoles() - protonHoles;
}


// Write information to output

std::ostream& operator<<(std::ostream& os, const G4ExitonConfiguration& ex) {
  os << " Exitons: protons " << ex.protonQuasiParticles << " holes " 
     << ex.protonHoles << "; neutrons " << ex.neutronQuasiParticles
     << " holes " << ex.neutronHoles;

  return os;
}
     

