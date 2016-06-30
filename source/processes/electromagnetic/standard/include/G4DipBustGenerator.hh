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
// $Id: G4DipBustGenerator.hh 96934 2016-05-18 09:10:41Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:  G4DipBustGenerator
//
// Author: Vladimir Grichine    
// Creation date: 17 May 2011
//
// Modifications: 
// 
//
// Bremsstrahlung Angular Distribution Generation 
// suggested the dipole approximation in the rest frame of electron 
// busted in the laboratory frame.
//
// -------------------------------------------------------------------
//

#ifndef G4DipBustGenerator_h
#define G4DipBustGenerator_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VEmAngularDistribution.hh"

class G4DipBustGenerator : public G4VEmAngularDistribution
{

public:

  explicit G4DipBustGenerator(const G4String& name = "");

  virtual ~G4DipBustGenerator();

  virtual G4ThreeVector& SampleDirection(const G4DynamicParticle* dp,
                                         G4double out_energy, G4int Z,
                                         const G4Material* mat = nullptr) final;

  G4double PolarAngle(const G4double initial_energy,
		      const G4double final_energy,
		      const G4int Z);

  virtual void PrintGeneratorInformation() const final;

private:

  // hide assignment operator 
  G4DipBustGenerator & operator=(const  G4DipBustGenerator &right) = delete;
  G4DipBustGenerator(const  G4DipBustGenerator&) = delete;

};

#endif

