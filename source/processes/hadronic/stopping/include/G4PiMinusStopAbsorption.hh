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
//      File name:     G4PiMinusStopAbsorption.hh 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 12 May 1998
//
//      Modifications: 
//      13 Sep 1998 - Changed DoAbsorption
//
// -------------------------------------------------------------------

#ifndef G4PIMINUSSTOPABSORPTION_HH
#define G4PIMINUSSTOPABSORPTION_HH

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4DynamicParticleVector.hh"
#include "G4PiMinusStopMaterial.hh"
#include "G4ThreeVector.hh"

class G4PiMinusStopAbsorption
{  

public:

  // Constructor
  G4PiMinusStopAbsorption(G4PiMinusStopMaterial* materialAlgo, const G4double Z, const G4double A);

  // Destructor
  ~G4PiMinusStopAbsorption();

  // Return final absorption products
  G4DynamicParticleVector* DoAbsorption();

  // Energy involved in the absorption process
  G4double Energy();

  // Return nucleus recoil momentum
  G4ThreeVector RecoilMomentum();

  // Number of protons in the absorption products
  G4int NProtons();

  // Number of neutrons in the absorption products
  G4int NNeutrons();

  void SetVerboseLevel(G4int level);

private:

  // Hide assignment operator as private 
  G4PiMinusStopAbsorption& operator=(const G4PiMinusStopAbsorption &right);

  // Copy constructor
  G4PiMinusStopAbsorption(const G4PiMinusStopAbsorption& );

  G4PiMinusStopMaterial* _materialAlgo; // owned pointer
  G4DynamicParticleVector* _absorptionProducts;

  G4double _nucleusA;
  G4double _nucleusZ;

  G4int _level;
};
 
#endif




















