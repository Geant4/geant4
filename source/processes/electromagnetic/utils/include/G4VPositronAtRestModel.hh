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
// File name:     G4VPositronAtRestModel
//
// Author:        V. Ivanchenko
// 
// Creation date: 13 May 2024
//
// Class Description: 
// Abstract base class for sampling of positron annihilation at rest
// There is a requirement to implementation of the virtual method
// SampleSecondaries(...):
//   in the list of dynamic particles two first should be the most
//   energetic gamma from the annihilation. The number of produced
//   particles is not limited.
//
// -------------------------------------------------------------------
//

#ifndef G4VPositronAtRestModel_h
#define G4VPositronAtRestModel_h 1

#include "globals.hh"
#include <vector>

class G4Material;
class G4DynamicParticle;

class G4VPositronAtRestModel
{
public:

  explicit G4VPositronAtRestModel(const G4String& name) : fName(name) {};

  virtual ~G4VPositronAtRestModel() = default;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>& secParticles,
				 G4double& localEnergyDeposit,
				 const G4Material*) const = 0;

  virtual void PrintGeneratorInformation() const = 0;

  const G4String& GetName() const { return fName; };

  G4VPositronAtRestModel& operator=
  (const  G4VPositronAtRestModel& right) = delete;
  G4VPositronAtRestModel(const G4VPositronAtRestModel&) = delete;

private:

  G4String fName;
};

#endif

