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
// $Id: G4DNARuddAngle.hh 66241 2012-12-13 18:34:42Z gunter $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:  G4DNARuddAngle
//
// Author:     Vladimir Ivantchenko
// 
// Creation date: 12 May 2014
//
// Modifications: 
//
// Class Description: 
//
// Delta-electron Angular Distribution Generation 
//
// -------------------------------------------------------------------
//

#ifndef G4DNARuddAngle_h
#define G4DNARuddAngle_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VEmAngularDistribution.hh"

class G4ParticleDefinition;

class G4DNARuddAngle : public G4VEmAngularDistribution
{

public:

  G4DNARuddAngle(const G4String& name = "");

  virtual ~G4DNARuddAngle();

  virtual G4ThreeVector& SampleDirection(const G4DynamicParticle* dp,
                                         G4double kinEnergyFinal,
                                         G4int Z,
                                         const G4Material* mat = 0);

  virtual G4ThreeVector& SampleDirectionForShell(
                                         const G4DynamicParticle* dp,
                                         G4double kinEnergyFinal,
                                         G4int Z, G4int shellIdx,
                                         const G4Material* mat = 0);

  void PrintGeneratorInformation() const;

private:

  // hide assignment operator 
  G4DNARuddAngle & operator=(const  G4DNARuddAngle &right);
  G4DNARuddAngle(const  G4DNARuddAngle&);

  const G4ParticleDefinition* fElectron;
};

#endif

