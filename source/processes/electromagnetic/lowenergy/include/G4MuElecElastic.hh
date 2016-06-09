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
// G4MuElecElastic.hh, 2011/08/29 A.Valentin, M. Raine
//
// Based on the following publications
//	    - Geant4 physics processes for microdosimetry simulation:
//	    very low energy electromagnetic models for electrons in Si,
//	    NIM B, vol. 288, pp. 66 - 73, 2012.
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef G4MuElecElastic_h
#define G4MuElecElastic_h 1

#include "G4VEmProcess.hh"
#include "G4Electron.hh"

// Available models
#include "G4MuElecElasticModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MuElecElastic : public G4VEmProcess

{
public: 

  G4MuElecElastic(const G4String& processName ="MuElecElastic",
		     G4ProcessType type = fElectromagnetic);

  virtual ~G4MuElecElastic();

  virtual G4bool IsApplicable(const G4ParticleDefinition&);
  
  virtual void PrintInfo();

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*);

private:
     
  G4bool       isInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
#endif
