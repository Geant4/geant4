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
// $Id: G4GammaConversion.hh 106717 2017-10-20 09:41:27Z gcosmo $
//
//
//------------------ G4GammaConversion physics process------------------------
//                   by Michel Maire, 24 May 1996
//
// 11-06-96, Added GetRandomAtom() method and new data member
//           for cumulative total cross section, by M.Maire
// 21-06-96, SetCuts inplementation, M.Maire
// 16-09-96, Dynamical array PartialSumSigma, M.Maire
// 14-01-97, crossection table + meanfreepath table.
//           PartialSumSigma removed, M.Maire
// 14-03-97, new physics scheme for geant4alpha, M.Maire
// 13-08-98, new methods SetBining() PrintInfo()
// 03-08-01, new methods Store/Retrieve PhysicsTable (mma)
// 06-08-01, BuildThePhysicsTable() called from constructor (mma)
// 19-09-01, come back to previous ProcessName: "conv"
// 20-09-01, DoIt: fminimalEnergy = 1*eV (mma)
// 01-10-01, come back to BuildPhysicsTable(const G4ParticleDefinition&)
// 13-08-04, suppress .icc file
//           public ComputeCrossSectionPerAtom() and ComputeMeanFreePath() (mma)
// 09-11-04, Remove Retrieve tables (V.Ivantchenko)
// 19-04-05, Redesign - use G4VEmProcess interface (V.Ivantchenko)
// 04-05-05, Make class to be default (V.Ivanchenko)
// 09-08-06, add SetModel(G4VEmModel*) (mma)
// 12-09-06, move SetModel(G4VEmModel*) in G4VEmProcess (mma)
// -----------------------------------------------------------------------------

// class description
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4GammaConversion_h
#define G4GammaConversion_h 1

#include "globals.hh"
#include "G4VEmProcess.hh"
#include "G4Gamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleDefinition;
class G4VEmModel;
class G4MaterialCutsCouple;
class G4DynamicParticle;

class G4GammaConversion : public G4VEmProcess

{
public:  // with description

  explicit G4GammaConversion(const G4String& processName ="conv",
			     G4ProcessType type = fElectromagnetic);

  virtual ~G4GammaConversion();

  // true for Gamma only.
  virtual G4bool IsApplicable(const G4ParticleDefinition&) final;

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
				    const G4Material*) override;

  // Print few lines of informations about the process: validity range,
  virtual void PrintInfo() override;

  // print documentation in html format
  virtual void ProcessDescription(std::ostream&) const override;

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*) override;

private:
     
  G4GammaConversion & operator=(const  G4GammaConversion &right) = delete;
  G4GammaConversion(const  G4GammaConversion&) = delete;

  G4bool  isInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
#endif
 
