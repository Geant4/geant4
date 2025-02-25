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
// File name:  G4RiGeAngularGenerator
//
// Authors:       Girardo Depaola & Ricardo Pacheco
// 
// Creation date: 29 October 2024
//
// -------------------------------------------------------------------
//

#ifndef G4RiGeAngularGenerator_h
#define G4RiGeAngularGenerator_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VEmAngularDistribution.hh"
#include "G4LorentzVector.hh"

class G4RiGeAngularGenerator : public G4VEmAngularDistribution
{

public:

  G4RiGeAngularGenerator();

  ~G4RiGeAngularGenerator() override = default;

  G4ThreeVector& SampleDirection(const G4DynamicParticle* dp,
                                 G4double gEnergy, G4int Z,
                                 const G4Material* mat = nullptr) override;

  G4LorentzVector Sample5DPairDirections(const G4DynamicParticle* dp,
                                 G4ThreeVector& dirElectron,
                                 G4ThreeVector& dirPositron,
                                 const G4double gEnergy, const G4double q2,
                                 const G4double gMomentum,
                                 G4double muFinalMomentum,
                                 G4double muFinalEnergy,
                                 const G4double* randNumbs,
                                 const G4double* W);

  void PhiRotation(G4ThreeVector& dir, G4double phi);

  G4LorentzVector eDP2(G4double x1, G4double x2, G4double x3, G4double x4, G4double x5);

  G4LorentzVector pDP2(G4double x3, const G4LorentzVector& x6);

  void PrintGeneratorInformation() const override;

  // hide assignment operator 
  G4RiGeAngularGenerator& operator=(const G4RiGeAngularGenerator& right) = delete;
  G4RiGeAngularGenerator(const G4RiGeAngularGenerator&) = delete;

private:

  G4double SampleCosTheta(G4double primKinEnergy, G4double gEnergy, G4double mass);

};

#endif

