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
// Geant4 header G4NuclearRadii
//
// Author V.Ivanchenko 27.05.2019
//
// Collection of parameterisations of nuclear radii selected from 
// different classes in cross section and model sub-libraries
//

#ifndef G4NuclearRadii_h
#define G4NuclearRadii_h 1

#include "globals.hh"

class G4Pow;
class G4ParticleDefinition;

class G4NuclearRadii 
{
public:

  // explicit radii for light nuclei
  static G4double ExplicitRadius(G4int Z, G4int A);

  // algorithm from diffuse-elastic parameterisation (V.Grichine)
  static G4double Radius(G4int Z, G4int A);

  // algorithm from e-A scattering data (V.Grichine)
  static G4double RadiusRMS(G4int Z, G4int A);

  // algorithm from Glauber-Gribov nucluear-nuclear model (V.Grichine)
  static G4double RadiusNNGG(G4int Z, G4int A);

  // algorithm of Edward Simpsons & Dousatsu Sakata
  static G4double RadiusECS(G4int Z, G4int A);

  // algorithm from Glauber-Gribov hadron-nuclear model (V.Grichine)
  static G4double RadiusHNGG(G4int A);

  // algorithm from Glauber-Gribov kaon-nuclear model (V.Grichine)
  static G4double RadiusKNGG(G4int A);

  // algorithm from nuclear de-excitation module
  static G4double RadiusND(G4int A);

  // algorithm from computation of Coulomb barrier in the nuclear 
  // de-excitation module
  static G4double RadiusCB(G4int Z, G4int A);

  // algorithm from computation of Coulomb barrier  
  static G4double ParticleRadius(const G4ParticleDefinition*);

  // algorithm for hadron-nucleon x-section
  static G4double CoulombFactor(
         const G4ParticleDefinition* theParticle, 
	 const G4ParticleDefinition* nucleon, 
	 G4double ekin);

  // algorithm for hadron-nucleus x-section
  static G4double CoulombFactor(
	 G4int Z, G4int A,
         const G4ParticleDefinition* theParticle, 
	 G4double ekin);

  static G4Pow* fG4pow;
  static const G4double r0[93];

};

#endif


