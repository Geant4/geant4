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

#ifndef TEST_HH
#define TEST_HH

#include "globals.hh"
#include <ostream>


class Test {

 public:
   Test(const G4String& physQuantity, 
        const G4String& category,
        G4double refValue, 
	G4double refRelError);
   ~Test();

   void SetSimulationValue(G4double value);

   void SetMetaData(G4double energy, 
                    const G4String& particle, 
                    const G4String& material) {

      primaryEnergy = energy;
      particleName = particle;
      materialName = material;
   }

   G4String GetPhysQuantity() {

      return quantity;
   }

   std::ostream& Print(std::ostream& os) const;

   G4bool Passed() { return passed; }

 private:
   const G4String quantity;
   const G4String unitCategory;

   const G4double referenceValue;
   const G4double referenceRelError;

   G4double computedValue;
   G4double diffPercent;
   G4int nmbSigmas;
   G4bool passed;

   G4double primaryEnergy;
   G4String particleName;
   G4String materialName;
};


inline std::ostream& operator<<(std::ostream& os, const Test& test) {

  return test.Print( os ); 
}

#endif
