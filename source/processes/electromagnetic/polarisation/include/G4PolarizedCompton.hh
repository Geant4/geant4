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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4PolarizedCompton
//
// Author:        Andreas Schaelicke
//                based on code by Michel Maire / Vladimir IVANTCHENKO
// Class description
//
// modified version respecting media and beam polarization
//     using the stokes formalism
//
// Creation date: 01.05.2005
//
// Modifications:
//
// 01-01-05, include polarization description (A.Stahl)
// 01-01-05, create asymmetry table and determine interactionlength (A.Stahl)
// 01-05-05, update handling of media polarization (A.Schalicke)
// 01-05-05, update polarized differential cross section (A.Schalicke)
// 26-07-06, cross section recalculated (P.Starovoitov)
// 09-08-06, make it work under current geant4 release (A.Schalicke)
// 11-06-07, add PostStepGetPhysicalInteractionLength (A.Schalicke)
//
// -----------------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PolarizedComptonScattering_h
#define PolarizedComptonScattering_h 1

#include "globals.hh"
#include "G4VEmProcess.hh"
#include "G4Gamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleDefinition;
class G4MaterialCutsCouple;
class G4DynamicParticle;

class G4PolarizedComptonModel;

class G4PolarizedCompton : public G4VEmProcess

{
public:  // with description

  explicit G4PolarizedCompton(const G4String& processName ="pol-compt",
		     G4ProcessType type = fElectromagnetic);

  virtual ~G4PolarizedCompton();

  // true for Gamma only.  
  virtual G4bool IsApplicable(const G4ParticleDefinition&) override;

  // Print few lines of informations about the process: validity range,
  virtual void PrintInfo() override;

  void SetModel(const G4String& name);

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*) override;

  // added for polarization treatment of polarized media:
  virtual void BuildPhysicsTable(const G4ParticleDefinition&) override;

  virtual G4double GetMeanFreePath(const G4Track& aTrack,     
				   G4double   previousStepSize,
				   G4ForceCondition* condition) override;

  virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition) override;

private:

  void CleanTable();

  void BuildAsymmetryTable(const G4ParticleDefinition& part);

  G4double ComputeAsymmetry(G4double energy,
			    const G4MaterialCutsCouple* couple,
			    const G4ParticleDefinition& particle,
			    G4double cut,
			    G4double & tAsymmetry);

  G4double ComputeSaturationFactor(const G4Track& aTrack);
  
  G4PolarizedCompton& operator=(const G4PolarizedCompton &right) = delete;
  G4PolarizedCompton(const G4PolarizedCompton& ) = delete;
     
  G4bool          buildAsymmetryTable;
  G4bool          useAsymmetryTable;

  G4bool          isInitialised;
  G4int           mType;

  // added for polarization treatment:
  G4PolarizedComptonModel* emModel;
  static G4PhysicsTable* theAsymmetryTable;  // table for crosssection assymmetry
  G4ThreeVector targetPolarization;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
#endif
 
