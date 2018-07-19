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
// $Id: G4DeltaAngle.hh 66241 2012-12-13 18:34:42Z gunter $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:  G4DeltaAngle
//
// Author:     Vladimir Ivantchenko
// 
// Creation date: 23 August 2013
//
// Modifications: 
//
// Class Description: 
//
// Delta-electron Angular Distribution Generation 
//
// -------------------------------------------------------------------
//

#ifndef G4DeltaAngle_h
#define G4DeltaAngle_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VEmAngularDistribution.hh"

class G4ParticleDefinition;

class G4DeltaAngle : public G4VEmAngularDistribution
{

public:

  explicit G4DeltaAngle(const G4String& name = "");

  virtual ~G4DeltaAngle();

  virtual G4ThreeVector& SampleDirection(const G4DynamicParticle* dp,
                                         G4double kinEnergyFinal, G4int Z,
                                         const G4Material* mat = nullptr) final;

  virtual G4ThreeVector& SampleDirectionForShell(
                                         const G4DynamicParticle* dp,
                                         G4double kinEnergyFinal,
                                         G4int Z, G4int shellIdx,
                                         const G4Material* mat = nullptr) final;

  virtual void PrintGeneratorInformation() const final;

private:

  // hide assignment operator 
  G4DeltaAngle & operator=(const  G4DeltaAngle &right) = delete;
  G4DeltaAngle(const  G4DeltaAngle&) = delete;

  const G4ParticleDefinition* fElectron;
  std::vector<G4double> prob;
  G4int                 nprob;
  G4int                 fShellIdx;
};

#endif

