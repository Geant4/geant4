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
// $Id: G4ICRU49NuclearStoppingModel.hh 103955 2017-05-04 11:29:54Z gcosmo $
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4ICRU49NuclearStoppingModel
//
// Author:        V.Ivanchenko 
//
// Creation date: 20.07.2009 
//
// Modifications:
//
//
// Class Description:
//
// Implementation of the ICRU'49 model of nuclear stopping 

// -------------------------------------------------------------------
//

#ifndef G4ICRU49NuclearStoppingModel_h
#define G4ICRU49NuclearStoppingModel_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4VEmModel.hh"
#include "G4Threading.hh"

class G4ParticleChangeForLoss;
class G4Pow;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ICRU49NuclearStoppingModel : public G4VEmModel
{

public:

  explicit G4ICRU49NuclearStoppingModel(const G4String& nam = 
					"ICRU49NucStopping");

  virtual ~G4ICRU49NuclearStoppingModel();

  virtual void Initialise(const G4ParticleDefinition*, 
			  const G4DataVector&) final;

  // main method to compute dEdx
  virtual G4double ComputeDEDXPerVolume(const G4Material*,
                                        const G4ParticleDefinition*,
                                        G4double kineticEnergy,
                                        G4double cutEnergy = DBL_MAX) final;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*, 
				 G4double, G4double) final;

private:

  void InitialiseArray();

  G4double NuclearStoppingPower(G4double kineticEnergy,
				G4double Z1, G4double Z2,
				G4double A1, G4double A2);

  //  hide assignment operator
  G4ICRU49NuclearStoppingModel & operator=
  (const  G4ICRU49NuclearStoppingModel &right) = delete;
  G4ICRU49NuclearStoppingModel(const  G4ICRU49NuclearStoppingModel&) = delete;

  G4Pow* g4calc;
  G4double theZieglerFactor;
  static G4double Z23[100];

#ifdef G4MULTITHREADED
  static G4Mutex ICRU49NuclearMutex;
#endif
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

