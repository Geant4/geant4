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
// GEANT4 Class file
//
//
// File name:     G4hQAOModel
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
//
// Creation date: 27 April 2004
//
// Modifications:
//
// Class Description: 
//
// Electronic stopping power for negative heavy partiles
// Further documentation available from http://www.ge.infn.it/geant4/lowE//

// -------------------------------------------------------------------
//

#ifndef G4hQAOModel_h
#define G4hQAOModel_h 1

#include "globals.hh"
#include "G4VhElectronicStoppingPower.hh"

class G4Material;
class G4Element;

class G4hQAOModel : public G4VhElectronicStoppingPower
{

public:

  G4hQAOModel();

  ~G4hQAOModel();

  G4bool HasMaterial(const G4Material*) {return true;};

  G4double StoppingPower(const G4Material* material,
                               G4double kineticEnergy);

  G4double ElectronicStoppingPower(G4double z,
                                   G4double kineticEnergy) const;

private:

  // get number of shell, energy and oscillator strenghts for material
  G4int GetNumberOfShell(G4int z) const;

  G4double GetShellEnergy(G4int z, G4int nbOfTheShell) const;
  G4double GetOscillatorEnergy(G4int z, G4int nbOfTheShell) const;
  G4double GetShellStrength(G4int z, G4int nbOfTheShell) const;
  G4double GetOccupationNumber(G4int z, G4int ShellNb) const;

  // calculate stopping number for L's term
  G4double GetL0(G4double normEnergy) const;
  // terms in Z^2
  G4double GetL1(G4double normEnergy) const;
  // terms in Z^3
  G4double GetL2(G4double normEnergy) const;
  // terms in Z^4
  
  // hide assignment operator
  G4hQAOModel & operator=(const  G4hQAOModel &right);
  G4hQAOModel(const  G4hQAOModel&);

  // Z of element at now avaliable for the model
  static const G4int materialAvailable[6];

  // number, energy and oscillator strenghts
  // for an harmonic oscillator model of material
  static const  G4int nbofShellForMaterial[6];
  static const  G4double alShellEnergy[3];
  static const  G4double alShellStrength[3];
  static const  G4double siShellEnergy[3];
  static const  G4double siShellStrength[3];
  static const  G4double cuShellEnergy[4];
  static const  G4double cuShellStrength[4];
  static const  G4double taShellEnergy[6];
  static const  G4double taShellStrength[6];
  static const  G4double auShellEnergy[6];
  static const  G4double auShellStrength[6];
  static const  G4double ptShellEnergy[6];
  static const  G4double ptShellStrength[6];


  //  variable for calculation of stopping number of L's term
  static const G4double L0[67][2];
  static const G4double L1[22][2];
  static const G4double L2[14][2];
  static const G4int nbOfElectronPerSubShell[1540];
  static const G4int fNumberOfShells[101];

  G4int numberOfMaterials;
  G4int sizeL0;
  G4int sizeL1;
  G4int sizeL2;

  const G4Material* currentMaterial;
  const G4Element*  currentElement;

  G4double theZieglerFactor;
  G4double thePlasmonFactor;
};

#endif
