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

#include "TestSeries.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "Test.hh"
#include "G4Material.hh"
#include <iomanip>


TestSeries::TestSeries(const G4String& quantity, 
                       const G4String& category, 
                       DetectorConstruction* det,
                       PrimaryGeneratorAction* gen) :
  currentTest(0),
  physQuantity(quantity),
  unitCategory(category),
  detector(det),
  particleGenerator(gen) {

}


TestSeries::~TestSeries() {

  std::vector<Test*>::iterator iter = tests.begin();
  std::vector<Test*>::iterator iter_end = tests.end();

  for(;iter != iter_end; iter++) delete (*iter);
  tests.clear();
} 


void TestSeries::CreateTestForRun(G4double refValue, G4double refRelError) {

  currentTest
         = new Test(physQuantity, unitCategory, refValue, refRelError);

  tests.push_back( currentTest );
}


void TestSeries::EndOfRunAction(G4double value) {

  if(currentTest != 0) {

     currentTest -> SetSimulationValue(value);

     G4double energy = particleGenerator -> 
       GetParticleGun()-> GetParticleEnergy();

     G4String particleName = particleGenerator -> GetParticleGun() ->
       GetParticleDefinition() -> GetParticleName();

     G4String materialName = detector -> GetAbsorMaterial() -> GetName();

     currentTest -> SetMetaData(energy,
                                particleName,
                                materialName);
  }

  // Test is reset
  currentTest = 0;
}


std::ostream& TestSeries::Print(std::ostream& os) const {

  os << G4endl;
  os << "==============================================================================="
     << G4endl;
  os << "                        TEST SUMMARY for  " << physQuantity
     << G4endl; 
  os << "  Ref.Value = Reference value, Err. = Relative error of reference" 
     << G4endl; 
  os << "  Sim.Value = Simulated value, Diff. = Relative difference Sim.-Ref. values" 
     << G4endl; 
  os << "  P = PASSED " 
     << G4endl; 
  os << "  F = FAILED " 
     << G4endl; 

  os << "-------------------------------------------------------------------------------"
     << G4endl;
  os << std::setw(10) << std::right << "Particle"
     << std::setw(15) << std::right << "Material"
     << std::setw(12) << std::right << "Energy"
     << std::setw(13) << std::right << "Ref.Value"
     << std::setw(6) << std::right << "Err."
     << std::setw(12) << std::right << "Sim.Value"
     << std::setw(7)  << std::right << "Diff."
     << std::setw(4)  << std::right << "P/F"
     << G4endl;
  os << std::setw(10) << std::right << ""
     << std::setw(15) << std::right << ""
     << std::setw(12) << std::right << "(MeV)"
     << std::setw(13) << std::right << ""
     << std::setw(6) << std::right << "(%)"
     << std::setw(12) << std::right << ""
     << std::setw(7)  << std::right << "(%)"
     << std::setw(4)  << std::right << ""
     << G4endl;
  os << "-------------------------------------------------------------------------------"
     << G4endl;

  std::vector<Test*>::const_iterator iter = tests.begin();
  std::vector<Test*>::const_iterator iter_end = tests.end();

  size_t passed = 0;
  size_t failed = 0;

  for(;iter != iter_end; iter++) { 

      os << *(*iter) << G4endl;

      if( (*iter) -> Passed() ) passed++;
      else failed ++;
  }

  os << "-------------------------------------------------------------------------------"
     << G4endl;
  os << "  TOTAL = " << passed + failed << ", FAILED = " << failed
     << ", PASSED = " << passed << G4endl;

  os << "==============================================================================="
     << G4endl;
  os << G4endl;

  return os;
}
