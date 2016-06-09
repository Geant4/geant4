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
// 20110922  M. Kelsey -- Replace print() with stream operator<<

#ifndef G4FISSION_CONFIGURATION_HH
#define G4FISSION_CONFIGURATION_HH

#include "globals.hh"
#include <iosfwd>

class G4FissionConfiguration {
public:
  G4FissionConfiguration() {}
  G4FissionConfiguration(G4double a, G4double z, G4double ez, 
			 G4double ek, G4double ep) 
    : afirst(a), zfirst(z), ezet(ez), ekin(ek), epot(ep) {}

  G4double afirst;
  G4double zfirst;
  G4double ezet;
  G4double ekin;
  G4double epot;
};        

std::ostream& operator<<(std::ostream& os, const G4FissionConfiguration& fis);

#endif // G4FISSION_CONFIGURATION_HH 

