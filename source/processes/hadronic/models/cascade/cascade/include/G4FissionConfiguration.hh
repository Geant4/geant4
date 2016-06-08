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
#ifndef G4FISSION_CONFIGURATION_HH
#define G4FISSION_CONFIGURATION_HH

#include "globals.hh"
#include "g4std/iostream"

class G4FissionConfiguration {

public:

  G4FissionConfiguration() {
  };

  G4FissionConfiguration(G4double a, 
			 G4double z, 
			 G4double ez, 
			 G4double ek, 
			 G4double ep) 
    : afirst(a), 
    zfirst(z), 
    ezet(ez), 
    ekin(ek), 
    epot(ep) {
};

  void print() {
    G4cout << " new configuration " << G4endl
	   << " a1 " << afirst << " z1 " << zfirst << " ez " << ezet <<
      " ekin " << ekin << " epot " << epot << G4endl;
  };

  G4double afirst;
  G4double zfirst;
  G4double ezet;
  G4double ekin;
  G4double epot;

};        

#endif // G4FISSION_CONFIGURATION_HH 

