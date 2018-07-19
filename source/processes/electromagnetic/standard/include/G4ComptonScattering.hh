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
// $Id: G4ComptonScattering.hh 106717 2017-10-20 09:41:27Z gcosmo $
//
//------------------ G4ComptonScattering physics process -----------------------
//                   by Michel Maire, April 1996
//
// 10-06-96, updated by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 06-01-97, crossection table + meanfreepath table, M.Maire
// 17-02-97, New Physics scheme
// 25-02-97, GetMeanFreePath() now is public function
// 12-03-97, new physics scheme again
// 13-08-98, new methods SetBining()  PrintInfo()
// 03-08-01, new methods Store/Retrieve PhysicsTable (mma)
// 06-08-01, BuildThePhysicsTable() called from constructor (mma)
// 19-09-01, come back to previous ProcessName "compt"
// 20-09-01, DoIt: fminimalEnergy = 1*eV (mma)
// 01-10-01, come back to BuildPhysicsTable(const G4ParticleDefinition&)
// 13-08-04, suppress icc file; make public ComputeCrossSectionPerAtom()   (mma)
// 09-11-04, Remove Retrieve tables (V.Ivantchenko)
// 15-03-05, Redesign - use G4VEmProcess interface (V.Ivantchenko)
// 04-05-05, Make class to be default (V.Ivanchenko)
// 09-09-06, modify SetModel(G4VEmModel*) (mma)
// 12-09-06, move SetModel(G4VEmModel*) in G4VEmProcess (mma)
// -----------------------------------------------------------------------------

// class description
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4ComptonScattering_h
#define G4ComptonScattering_h 1

#include "globals.hh"
#include "G4VEmProcess.hh"
#include "G4Gamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleDefinition;
class G4VEmModel;
class G4MaterialCutsCouple;
class G4DynamicParticle;

class G4ComptonScattering : public G4VEmProcess

{
public:  // with description

  explicit G4ComptonScattering(const G4String& processName ="compt",
			       G4ProcessType type = fElectromagnetic);

  virtual ~G4ComptonScattering();

  // true for Gamma only.  
  virtual G4bool IsApplicable(const G4ParticleDefinition&) final;
  
  // Print few lines of informations about the process: validity range,
  virtual void PrintInfo() override;
 
  // print description in html
  virtual void ProcessDescription(std::ostream&) const override;

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*) override;

private:
     
  G4ComptonScattering & operator=(const  G4ComptonScattering &right) = delete;
  G4ComptonScattering(const  G4ComptonScattering&) = delete;

  G4bool       isInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
#endif
 
