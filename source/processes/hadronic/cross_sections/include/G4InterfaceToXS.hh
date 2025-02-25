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
// GEANT4 Class header file
//
// File name:  G4InterfaceToXS
//
// Author V.Ivantchenko, 12 July 2024
//
//
// Class description:
// Interface to G4PARTICLEXS cross sections for the de_excitation module
//

#ifndef G4InterfaceToXS_h
#define G4InterfaceToXS_h 1

#include "globals.hh"

class G4GammaNuclearXS;
class G4NeutronInelasticXS;
class G4ParticleInelasticXS;
class G4ParticleDefinition;

class G4InterfaceToXS
{
public:

  G4InterfaceToXS(const G4ParticleDefinition*, G4int index);

  ~G4InterfaceToXS() = default;

  void Initialise();
  
  G4double GetElementCrossSection(const G4double ekin, const G4int Z); 

  G4double GetIsoCrossSection(const G4double ekin, const G4int Z, const G4int A); 

  G4InterfaceToXS(const G4InterfaceToXS& right) = delete;
  const G4InterfaceToXS& operator = (const G4InterfaceToXS& right) = delete;

private:

  G4int index;
  const G4ParticleDefinition* fParticle;
  G4GammaNuclearXS* fGammaNuclear{nullptr};
  G4NeutronInelasticXS* fNeutronNuclear{nullptr};
  G4ParticleInelasticXS* fParticleNuclear{nullptr};
};

#endif
