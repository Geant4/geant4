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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:  G4ModifiedMephi
//
// Author:        V. Ivanchenko
// 
// Creation date: 27 October 2020
//
// Modifications: 
//
// Class Description: 
//
// Bremsstrahlung Angular Distribution Generation 
// A.G. Bogdanov et al., IEEE Trans. Nuc. Sci., Vol.53, No.2, 2006
//
// -------------------------------------------------------------------
//

#ifndef G4ModifiedMephi_h
#define G4ModifiedMephi_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VEmAngularDistribution.hh"

class G4ModifiedMephi : public G4VEmAngularDistribution
{

public:

  explicit G4ModifiedMephi(const G4String& name = "");

  virtual ~G4ModifiedMephi();

  virtual G4ThreeVector& SampleDirection(const G4DynamicParticle* dp,
                                         G4double gEnergy, G4int Z,
                                         const G4Material* mat = nullptr) final;

  virtual void SamplePairDirections(const G4DynamicParticle* dp,
				    G4double elecKinEnergy,
				    G4double posiKinEnergy,
				    G4ThreeVector& dirElectron,
				    G4ThreeVector& dirPositron,
				    G4int Z = 0,
				    const G4Material* mat = nullptr) final;

  virtual void PrintGeneratorInformation() const final;

  // hide assignment operator 
  G4ModifiedMephi & operator=(const  G4ModifiedMephi &right) = delete;
  G4ModifiedMephi(const  G4ModifiedMephi&) = delete;

private:

  G4double SampleCosTheta(G4double primKinEnergy, G4double gEnergy,
                          G4double mass);

};

#endif

