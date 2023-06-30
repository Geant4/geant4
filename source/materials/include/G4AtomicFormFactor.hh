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

//---------------------------------------------------------------------------
//
// ClassName:   G4AtomicFormFactor
//
// Description: Contains the function for the evaluation of the atomic form
//              factor. The tabulated data are available on IUCr website
//
// Class description:
//
// XXX
//
// 21-04-16, created by E.Bagli

#ifndef G4ATOMICFORMFACTOR_HH
#define G4ATOMICFORMFACTOR_HH 1

#include "G4Exp.hh"
#include "globals.hh"

#include <map>
#include <vector>

class G4AtomicFormFactor
{
 public:
  static G4AtomicFormFactor* GetManager()
  {
    if (s_G4AtomicFormFactorManager == nullptr) {
      s_G4AtomicFormFactorManager = new G4AtomicFormFactor();
    }
    return s_G4AtomicFormFactorManager;
  }

  G4double Get(G4double kScatteringVector, G4int Z, G4int charge = 0)
  {
    if (loadedIndex != GetIndex(Z, charge)) {
      LoadCoefficiencts(GetIndex(Z, charge));
    }
    G4double result = 0.;
    G4double kVecOn4PiSquared =
      (kScatteringVector / 1.e-7 / 3.1415926536) * 0.125;  // (k/(4pi))/   angstrom
    kVecOn4PiSquared = kVecOn4PiSquared * kVecOn4PiSquared;  // (k/(4pi))^2

    for (unsigned int i0 = 0; i0 < 4; i0++) {
      result += theCoefficients[i0 * 2] * G4Exp(-theCoefficients[i0 * 2 + 1] * kVecOn4PiSquared);
    }
    result += theCoefficients[8];
    return result;
  }

 protected:
  G4AtomicFormFactor(); 
  ~G4AtomicFormFactor() = default;

 private:
  void InsertCoefficients(G4int index, const std::vector<G4double>& aDoubleVec)
  {
    theCoefficientsMap.insert(std::pair<G4int, std::vector<G4double>>(index, aDoubleVec));
  }

  // LoadCoefficiencts() method allows the evaluation of the atomic form
  // factor coefficients and the storage in theCoefficients.
  // If theCoefficients are already correct, no need to get new ones
  // Reference: International Tables for Crystallography (2006).
  // Vol. C, ch. 6.1, pp. 554-595
  // doi: 10.1107/97809553602060000600
  // Chapter 6.1. Intensity of diffracted intensities
  // IUCr Eq. 6.1.1.15, Coefficients Table 6.1.1.4
  void LoadCoefficiencts(G4int index)
  {
    loadedIndex = index;
    for (unsigned int i0 = 0; i0 < 9; i0++) {
      theCoefficients[i0] = theCoefficientsMap[index][i0];
    }
  }

  // Get() function gives back the Atomic Form Factor of the Z material
  inline G4int GetIndex(G4int Z, G4int charge = 0) { return Z * 100 + charge; }

 private:
  inline static G4AtomicFormFactor* s_G4AtomicFormFactorManager = nullptr;

  // theCoefficientsMap stores the coefficients for the form factor
  // calculations. It can be loaded only by LoadCoefficiencts()
  // and accessed by theCoefficients[].
  std::map<G4int, std::vector<G4double>> theCoefficientsMap;
  G4double theCoefficients[9];
  G4int loadedIndex;
};
#endif
