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
// $Id:  G4RayleighAngularGenerator.hh,v 1.1 2012-12-31 18:34:15 antoni Exp $
// GEANT4 tag $Name: not supported by svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:  G4RayleighAngularGenerator
//
// Author:     Ivantchenko  (antoni@cern.ch, vnivanch@cern.ch)
//             modified fit formulas from Dermott E. Cullen, 
//             Nucl. Instrum. Meth. Phys. Res. B v.101, (4),499-510. 
//            
// 
// Creation date: 31 May 2012
//
// Modifications: 
//
//
// Class Description: 
//
// Class for Rayleigh Scattering angle sampling
//
// -------------------------------------------------------------------
//

#ifndef G4RayleighAngularGenerator_h
#define G4RayleighAngularGenerator_h 1

#include "G4VEmAngularDistribution.hh"
#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"

class G4RayleighAngularGenerator : public G4VEmAngularDistribution
{
public:

  G4RayleighAngularGenerator();

  virtual ~G4RayleighAngularGenerator();

  virtual G4ThreeVector& SampleDirection(const G4DynamicParticle* dp,
					 G4double out_energy,
					 G4int Z,
					 const G4Material* mat = 0);
private:

  // hide assignment operator 
  G4RayleighAngularGenerator & operator=(const  G4RayleighAngularGenerator &right);
  G4RayleighAngularGenerator(const  G4RayleighAngularGenerator&);

  // static data
  static const G4double PP0[101];
  static const G4double PP1[101];
  static const G4double PP2[101];
  static const G4double PP3[101];
  static const G4double PP4[101];
  static const G4double PP5[101];
  static const G4double PP6[101];
  static const G4double PP7[101];
  static const G4double PP8[101];

  G4double fFactor;
  
};


#endif

