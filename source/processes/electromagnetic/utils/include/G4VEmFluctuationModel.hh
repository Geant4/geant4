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
// $Id: G4VEmFluctuationModel.hh 106208 2017-09-20 01:53:57Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4VEmFluctuationModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.01.2002
//
// Modifications:
//
// 28-12-02 add method Dispersion (V.Ivanchenko)
// 07-02-03 change signature (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 16-10-03 Changed interface to Initialisation (V.Ivanchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 25-07-05 Move constructor and destructor to the body (V.Ivanchenko)
//
//
// Class Description: 
//
// Abstract class for interface to simualtion of energy loss fluctuations

// -------------------------------------------------------------------
//

#ifndef G4VEmFluctuationModel_h
#define G4VEmFluctuationModel_h 1


#include "globals.hh"
#include <CLHEP/Random/RandomEngine.h>

class G4ParticleDefinition;
class G4DynamicParticle;
class G4Material;
class G4MaterialCutsCouple;
class G4LossTableManager;

class G4VEmFluctuationModel 
{

public:

  explicit G4VEmFluctuationModel(const G4String& nam);

  virtual ~G4VEmFluctuationModel();

  //------------------------------------------------------------------------
  // Virtual methods to be implemented for the concrete model
  //------------------------------------------------------------------------

  virtual G4double SampleFluctuations(const G4MaterialCutsCouple*,
				      const G4DynamicParticle*,
				      G4double tmax,
				      G4double length,
				      G4double meanLoss) = 0;

  virtual G4double Dispersion(const G4Material*,
                              const G4DynamicParticle*,
			      G4double tmax,
			      G4double length) = 0;

  //------------------------------------------------------------------------
  // Methods with standard implementation; may be overwritten if needed 
  //------------------------------------------------------------------------

  virtual void InitialiseMe(const G4ParticleDefinition*);

  virtual void SetParticleAndCharge(const G4ParticleDefinition*, G4double q2);

  //------------------------------------------------------------------------
  // Generic methods common to all models
  //------------------------------------------------------------------------

  inline const G4String& GetName() const;

private:

  // hide assignment operator
  G4VEmFluctuationModel & 
    operator=(const  G4VEmFluctuationModel &right) = delete;
  G4VEmFluctuationModel(const  G4VEmFluctuationModel&) = delete;

  const G4String      name;
  G4LossTableManager* fManager;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4String& G4VEmFluctuationModel::GetName() const 
{
  return name;
}

#endif

