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
#include "G4FissionConfiguration.hh"

#include "g4std/vector"

class G4FissionStore {

public:

  G4FissionStore();

  void addConfig(G4double a, 
		 G4double z, 
		 G4double ez, 
		 G4double ek, 
		 G4double ev) {
    G4FissionConfiguration config(a, z, ez, ek, ev);
    configurations.push_back(config);
    // config.print();
  };

  G4int size() const { 
    return configurations.size(); 
  };

  G4FissionConfiguration generateConfiguration(G4double amax, 
					       G4double rand) const;

private:
  G4int verboseLevel;

  G4std::vector<G4FissionConfiguration> configurations;

};

#endif // G4FISSION_STORE_HH 



