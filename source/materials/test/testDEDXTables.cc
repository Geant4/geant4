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
// Anton.Lechner@cern.ch
//
// Environment variable TESTTARGET must be specified (see GNUmakefile)
//

#include "G4MaterialStoppingICRU73.hh"
#include "G4SimpleMaterialStoppingICRU73.hh"
#include "G4IronStoppingICRU73.hh"
#include <sstream>


G4bool Compare(G4PhysicsVector* physicsVector,
               G4double refDEDXValue,
               G4double kinEnergyPerNucleon, 
               G4int atomicNumberIon,
               const G4String& matIdentifier) {

  G4double tableDEDXValue = -1.0;

  if(physicsVector != 0) {
     size_t nmbBins = physicsVector -> GetVectorLength();
     G4double lowerEnergyEdge = physicsVector -> GetLowEdgeEnergy(0);
     G4double upperEnergyEdge = physicsVector -> GetLowEdgeEnergy(nmbBins-1);

     if(kinEnergyPerNucleon <= upperEnergyEdge && 
        kinEnergyPerNucleon >= lowerEnergyEdge) {
  	G4bool b;
 
        tableDEDXValue = physicsVector -> GetValue(kinEnergyPerNucleon, b);
     }
  }

  //  G4cout << *physicsVector << G4endl;

  G4cout << " E/A1(MeV) = "   << kinEnergyPerNucleon / MeV
         << ", Z(Ion) = "     << atomicNumberIon
         << ", Mat = "        << matIdentifier
         << " -- dE/dx(MeV*cm2/mg) = "    
         << tableDEDXValue / (MeV * cm2 / (0.001 * g))
	 << ", Ref(MeV*cm2/mg) = " 
         << refDEDXValue / (MeV * cm2 / (0.001 * g))
         << ":";

  G4double eps = 0.0001 * (MeV * cm2 / (0.001 * g));

  if( refDEDXValue < (tableDEDXValue - eps) ||
      refDEDXValue > (tableDEDXValue + eps) ) {

     G4cout << " Test failed" << G4endl;
     return false;
  } 
  else G4cout << " Test passed" << G4endl; 

  return true;
}


G4bool CompareDEDXValues(G4VIonDEDXTable* table,
                         G4double refDEDXValue,
                         G4double kinEnergyPerNucleon, 
                         G4int atomicNumberIon,
                         G4int atomicNumberMat) {

  G4PhysicsVector* physicsVector = 
                            table -> GetPhysicsVector(atomicNumberIon,
                                                      atomicNumberMat);

  std::stringstream ss;
  ss << atomicNumberMat;

  return Compare(physicsVector, refDEDXValue, kinEnergyPerNucleon, 
                 atomicNumberIon, ss.str());
}


G4bool CompareDEDXValues(G4VIonDEDXTable* table,
                         G4double refDEDXValue,
                         G4double kinEnergyPerNucleon, 
                         G4int atomicNumberIon,
                         const G4String& matName) {

  G4PhysicsVector* physicsVector = 
                            table -> GetPhysicsVector(atomicNumberIon,
                                                      matName);
  return Compare(physicsVector, refDEDXValue, kinEnergyPerNucleon, 
                 atomicNumberIon, matName);
}


G4bool TestApplicability(G4VIonDEDXTable* table,
                         G4bool shouldBeApplicable,
                         G4int atomicNumberIon,
                         const G4String& matName) {

  G4bool isApplicable = table -> IsApplicable(atomicNumberIon, matName);

  G4cout << " Z(Ion) = " << atomicNumberIon
         << ", Mat = "   << matName
         << " -- ";

  if(shouldBeApplicable) {
     G4cout << " Should be applicable: "; 
     if(isApplicable) G4cout << "Test passed." << G4endl;
     else { G4cout << "Test failed." << G4endl; return false; }
  }
  else {
     G4cout << " Should not be applicable: ";   
     if(!isApplicable) G4cout << "Test passed." << G4endl;
     else { G4cout << "Test failed." << G4endl; return false; }
  }
  return true;
}


int main() {  

  G4double unitDEDX = (MeV * cm2 / (0.001 * g));

  G4SimpleMaterialStoppingICRU73* icru73elem = 
                                         new G4SimpleMaterialStoppingICRU73();

  G4MaterialStoppingICRU73* icru73comp = new G4MaterialStoppingICRU73();

  G4IronStoppingICRU73* icru73iron = new G4IronStoppingICRU73();

  G4cout << "### Test of DEDX tables (ICRU 73) ###" << G4endl; 

  G4cout << "#### I. Testing Function GetPhysicsVector() ###" << G4endl;

  G4cout << G4endl << "### A. G4SimpleMaterialStoppingICRU73:" << G4endl;
    
  // Reference values taken directly from ICRU 73 tables
  CompareDEDXValues(icru73elem, 8.31  * unitDEDX, 0.025, 3, 1 );
  CompareDEDXValues(icru73elem, 3.312 * unitDEDX, 1.5,   8, 50);
  CompareDEDXValues(icru73elem, 10.99 * unitDEDX, 10., 18, 10);
  CompareDEDXValues(icru73elem, 4.412 * unitDEDX, 10., 18, 92);

  CompareDEDXValues(icru73elem, 8.31  * unitDEDX, 0.025, 3, "G4_H");
  CompareDEDXValues(icru73elem, 3.312 * unitDEDX, 1.5,   8, "G4_Sn");
  CompareDEDXValues(icru73elem, 10.99 * unitDEDX, 10., 18, "G4_Ne");
  CompareDEDXValues(icru73elem, 4.412 * unitDEDX, 10., 18, "G4_U");

  G4cout << G4endl << "### B. G4MaterialStoppingICRU73:" << G4endl;
    
  // Reference values taken directly from ICRU 73 tables
  CompareDEDXValues(icru73comp, 2.748  * unitDEDX, 0.025, 3, "G4_A-150_TISSUE");
  CompareDEDXValues(icru73comp, 1.909  * unitDEDX, 0.025, 3, "G4_CARBON_DIOXIDE");
  CompareDEDXValues(icru73comp, 4.225  * unitDEDX, 10.0 , 10, "G4_MYLAR");

  // Reference values taken directly from revised ICRU 73 tables
  CompareDEDXValues(icru73comp, 0.19364 * unitDEDX, 25.0, 3, "G4_WATER");
  CompareDEDXValues(icru73comp, 3.6037 * unitDEDX, 0.025, 6, "G4_WATER");
  CompareDEDXValues(icru73comp, 1.6302 * unitDEDX, 10.0, 6, "G4_WATER");
  CompareDEDXValues(icru73comp, 0.079682 * unitDEDX, 1000.0, 6, "G4_WATER");
  CompareDEDXValues(icru73comp, 1.2296 * unitDEDX, 80.0, 12, "G4_WATER");
  CompareDEDXValues(icru73comp, 7.1917 * unitDEDX, 15.0, 15, "G4_WATER");
  CompareDEDXValues(icru73comp, 4.7505 * unitDEDX, 30.0, 16, "G4_WATER");
  CompareDEDXValues(icru73comp, 0.76464 * unitDEDX, 800.0, 18, "G4_WATER");
  CompareDEDXValues(icru73comp, 12.836 * unitDEDX, 0.09, 18, "G4_WATER");
  CompareDEDXValues(icru73comp, 0.39022 * unitDEDX, 250.0, 10, "G4_WATER");

  G4cout << G4endl << "### C. G4IronStoppingICRU73:" << G4endl;

  // Reference values taken directly from ICRU 73 tables
  CompareDEDXValues(icru73iron, 27.6  * unitDEDX, 0.025, 26, "G4_H");
  CompareDEDXValues(icru73iron, 30.72  * unitDEDX, 1.5, 26, "G4_Al");
  CompareDEDXValues(icru73iron, 13.75  * unitDEDX, 10.0 , 26, "G4_Kr");

  G4cout << G4endl << "#### II. Testing Function IsApplicable() ###" << G4endl;

  G4cout << G4endl << "### A. G4SimpleMaterialStoppingICRU73:" << G4endl;
    
  TestApplicability(icru73elem, true, 3, "G4_H");
  TestApplicability(icru73elem, true, 8, "G4_Sn");
  TestApplicability(icru73elem, true, 18, "G4_Ne");
  TestApplicability(icru73elem, true, 18, "G4_U");
  TestApplicability(icru73elem, false, 19, "G4_H");
  TestApplicability(icru73elem, false, 2, "G4_Sn");
  TestApplicability(icru73elem, false, 3, "G4_Cr");

  G4cout << G4endl << "### B. G4MaterialStoppingICRU73:" << G4endl;
    
  TestApplicability(icru73comp, true, 3, "G4_MYLAR");
  TestApplicability(icru73comp, true, 8, "G4_POLYETHYLENE");
  TestApplicability(icru73comp, true, 18, "G4_CARBON_DIOXIDE");
  TestApplicability(icru73comp, true, 18, "G4_A-150_TISSUE");
  TestApplicability(icru73comp, false, 19, "G4_MYLAR");
  TestApplicability(icru73comp, false, 2, "G4_MYLAR");
  TestApplicability(icru73comp, false, 3, "G4_Cr");

  G4cout << G4endl << "### B. G4IronStoppingICRU73:" << G4endl;
    
  TestApplicability(icru73iron, true, 26, "G4_H");
  TestApplicability(icru73iron, true, 26, "G4_Ne");
  TestApplicability(icru73iron, true, 26, "G4_Kr");
  TestApplicability(icru73iron, false, 19, "G4_H");
  TestApplicability(icru73iron, false, 26, "G4_Sn");
  TestApplicability(icru73iron, false, 3, "G4_Cr");

  G4cout << G4endl << "################# DONE #################" << G4endl;

  delete icru73elem;
  delete icru73comp;
  delete icru73iron;

  return 0;
} 

