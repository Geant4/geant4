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
// $Id: G4PhotoElectricEffect.hh 106717 2017-10-20 09:41:27Z gcosmo $
//
//------------------ G4PhotoElectricEffect physics process ------------------
//                   by Michel Maire, 24 May 1996
//
// 12-06-96, Added SelectRandomAtom() method and new data member
//           for cumulative total cross section, by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 17-09-96, Dynamic array PartialSumSigma
//           split ComputeBindingEnergy(), M.Maire
// 08-01-97, crossection table + meanfreepath table, M.Maire
// 13-03-97, adapted for the new physics scheme, M.Maire
// 13-08-98, new methods SetBining() PrintInfo()
// 17-11-98, use table of atomic shells in PostStepDoIt, mma
// 06-01-99, Sandia crossSection below 50 keV, V.Grichine mma
// 03-08-01, new methods Store/Retrieve PhysicsTable (mma)
// 06-08-01, BuildThePhysicsTable() called from constructor (mma)
// 19-09-01, come back to previous process name "phot"
// 20-09-01, DoIt: fminimalEnergy = 1*eV (mma)
// 01-10-01, come back to BuildPhysicsTable(const G4ParticleDefinition&)
// 10-01-02, moved few function from icc to cc
// 17-04-02, Keep only Sandia crossSections. Remove BuildPhysicsTables.
//           Simplify public interface (mma)
// 29-04-02, Generate theta angle of the photoelectron from Sauter-Gavrila
//           distribution (mma)
// 13-08-04, suppress icc file; make public ComputeCrossSectionPerAtom() (mma)
// 21-04-05, Redesign - use G4VEmProcess interface (V.Ivanchenko)
// 02-05-05, move ParticleChange actions in model (mma)
// 04-05-05, Make class to be default (V.Ivanchenko)
// 09-08-06, add SetModel(G4VEmModel*) (mma)
// 12-09-06, move SetModel(G4VEmModel*) in G4VEmProcess (mma)
// -----------------------------------------------------------------------------

// class description
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4PhotoElectricEffect_h
#define G4PhotoElectricEffect_h 1

#include "globals.hh"
#include "G4VEmProcess.hh"
#include "G4Gamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleDefinition;
class G4VEmModel;
class G4MaterialCutsCouple;
class G4DynamicParticle;

class G4PhotoElectricEffect : public G4VEmProcess

{
public:  // with description

  explicit G4PhotoElectricEffect(const G4String& processName ="phot",
				 G4ProcessType type = fElectromagnetic);

  virtual ~G4PhotoElectricEffect();

  // true for Gamma only.
  virtual G4bool IsApplicable(const G4ParticleDefinition&) final;

  // Print few lines of informations about the process: validity range,
  virtual void PrintInfo() override;

  // print documentation in html format
  virtual void ProcessDescription(std::ostream&) const override;

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*) override;

private:

  G4bool  isInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

