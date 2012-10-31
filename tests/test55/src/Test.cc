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

#include <iomanip>

#include "Test.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

Test::Test(const G4String& physQuantity, 
           const G4String& category, 
           G4double refValue, 
           G4double refRelError) :
   quantity(physQuantity), 
   unitCategory(category),
   referenceValue(refValue),
   referenceRelError(refRelError),
   computedValue(0),
   diffPercent(0),
   nmbSigmas(0),
   passed(false),
   primaryEnergy(0),
   particleName(""),
   materialName("") {

}


Test::~Test() {

}


void Test::SetSimulationValue(G4double value) {

  computedValue = value;

  diffPercent = (1.0 - computedValue/referenceValue) * 100.0;

  nmbSigmas = 0;

  G4double sigma = referenceValue * referenceRelError;
  G4double distance = sigma;

  if(computedValue <= referenceValue) {

    nmbSigmas--;

    while( computedValue < referenceValue - distance ) {

      nmbSigmas--;
      distance += sigma; 
    }
  }
  else {

    nmbSigmas++;

    while( computedValue > referenceValue + distance ) {

      nmbSigmas++;
      distance += sigma; 
    }
  }

  if(nmbSigmas == 1 || nmbSigmas == -1) passed = true;
}


std::ostream& Test::Print(std::ostream& os) const {

  std::ios::fmtflags mode = os.flags();

  os.setf(std::ios::fixed,std::ios::floatfield);

  G4long prec = os.precision(1);

  std::ostringstream ossRef;
  ossRef << G4BestUnit(referenceValue, unitCategory);

  std::ostringstream ossSim;
  ossSim << G4BestUnit(computedValue, unitCategory);

  G4String material = materialName;
  if(materialName.length() > 13) 
      material = materialName.substr(0,10) + "...";

  os << std::setw(10) << std::right << particleName
     << std::setw(15) << std::right << material
     << std::setw(12) << std::right << primaryEnergy/MeV
     << std::setw(14) << std::right << ossRef.str()
     << std::setw(5) << std::right << referenceRelError * 100.0 
     << std::setw(13) << std::right << ossSim.str()
    //     << std::setw(8)  << std::right << nmbSigmas;
     << std::setw(6)  << std::right << diffPercent;

  if(passed) os << std::setw(4)  << std::right << "P";
  else os << std::setw(4)  << std::right << "F";

  os.precision(prec);    
  os.setf(mode,std::ios::floatfield);

  return os;
}
