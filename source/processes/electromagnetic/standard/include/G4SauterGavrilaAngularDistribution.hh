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
// File name:  G4SauterGavrilaAngularDistribution
//
// Author:     Vladimir Ivanchenko using Michel Maire algorithm
//             developed for Geant3
// 
// Creation date: 23 July 2012
//
// Class Description: 
//
// Photoelectric Angular Distribution Generation 
// The Sauter-Gavrila distribution for the K-shell is used.
//
// -------------------------------------------------------------------
//

#ifndef G4SauterGavrilaAngularDistribution_h
#define G4SauterGavrilaAngularDistribution_h 1

#include "G4VEmAngularDistribution.hh"

class G4SauterGavrilaAngularDistribution : public G4VEmAngularDistribution
{

public:

  explicit G4SauterGavrilaAngularDistribution();

  ~G4SauterGavrilaAngularDistribution() override;

  G4ThreeVector& SampleDirection(const G4DynamicParticle* dp,
                                 G4double e = 0.0,
                                 G4int shellId = 0,
                                 const G4Material* mat = nullptr) final;

  void PrintGeneratorInformation() const override;

  // hide assignment operator 
  G4SauterGavrilaAngularDistribution & operator=
  (const  G4SauterGavrilaAngularDistribution &right) = delete;
  G4SauterGavrilaAngularDistribution(const  G4SauterGavrilaAngularDistribution&) = delete;

};

#endif

