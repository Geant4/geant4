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
// $Id: G4ePolarizedIonisation.hh 68046 2013-03-13 14:31:38Z gcosmo $
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4ePolarizedIonisation
//
// Author:        A.Schaelicke on base of Vladimir Ivanchenko code
//
// Creation date: 10.11.2005
//
// Modifications:
//
// 10-11-05, include polarization description (A.Schaelicke)
// , create asymmetry table and determine interactionlength 
// , update polarized differential cross section 
//
// 20-08-06, modified interface (A.Schaelicke)
// 11-06-07, add PostStepGetPhysicalInteractionLength (A.Schalicke)
//
// Class Description:
//
// polarized version of G4eIonisation
// ----------------------------------------------------------------------------
// -------------------------------------------------------------------
//

#ifndef G4ePolarizedIonisation_h
#define G4ePolarizedIonisation_h 1

#include "G4VEnergyLossProcess.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

class G4Material;
class G4ParticleDefinition;
class G4VEmFluctuationModel;
class G4PolarizedMollerBhabhaModel;

class G4ePolarizedIonisation : public G4VEnergyLossProcess
{

public:

  G4ePolarizedIonisation(const G4String& name = "pol-eIoni");

  virtual ~G4ePolarizedIonisation();

  G4bool IsApplicable(const G4ParticleDefinition& p);

  // Print out of the class parameters
  virtual void PrintInfo();

protected:

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
					   const G4ParticleDefinition*);

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
                                    const G4Material*, G4double cut);


  // for polarization
  G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            );

  G4double GetMeanFreePath(const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition);

  virtual void BuildPhysicsTable(const G4ParticleDefinition&);

  G4double ComputeAsymmetry(G4double energy,
			    const G4MaterialCutsCouple* couple,
			    const G4ParticleDefinition& particle,
			    G4double cut,
			    G4double &tasm);

protected:

  const G4ParticleDefinition* DefineBaseParticle(const G4ParticleDefinition* p);

private:
  G4ePolarizedIonisation & operator=(const G4ePolarizedIonisation &right);
  G4ePolarizedIonisation(const G4ePolarizedIonisation&);

  G4ParticleDefinition* theElectron;
  G4VEmFluctuationModel* flucModel;
  G4PolarizedMollerBhabhaModel* emModel;

  G4bool isElectron;
  G4bool isInitialised;

  // for polarization:
  G4ThreeVector theTargetPolarization;


  
  G4PhysicsTable* theAsymmetryTable;          // table for cross section assym.
  G4PhysicsTable* theTransverseAsymmetryTable; // table for transverse cross section assym.
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4ePolarizedIonisation::MinPrimaryEnergy(const G4ParticleDefinition*,
                                                   const G4Material*,
                                                         G4double cut)
{
  G4double x = cut;
  if(isElectron) x += cut;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4ePolarizedIonisation::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Electron::Electron() || &p == G4Positron::Positron());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
