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
// $Id: G4PolarizedGammaConversion.hh 96114 2016-03-16 18:51:33Z gcosmo $
//
// -------------------------------------------------------------------
//
//
// File name:     G4PolarizedGammaConversion
//
// Author:        Karim Laihem based on code by Michel Maire
//
// Creation date: 01.05.2005
//
// Modifications:
// 21-08-06 Modified to work in g4.8.1 framework (A.Schaelicke)
//
// Class Description:
//
// polarized version of G4GammaConversion

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4PolarizedGammaConversion_h
#define G4PolarizedGammaConversion_h 1

#include "globals.hh"
#include "G4VEmProcess.hh"
#include "G4Gamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleDefinition;
class G4VEmModel;
class G4MaterialCutsCouple;
class G4DynamicParticle;

class G4PolarizedGammaConversion : public G4VEmProcess

{
public:  // with description

  explicit G4PolarizedGammaConversion(const G4String& processName ="pol-conv",
			     G4ProcessType type = fElectromagnetic);

  virtual ~G4PolarizedGammaConversion();

  // true for Gamma only.
  G4bool IsApplicable(const G4ParticleDefinition&) override;

  // Print few lines of informations about the process: validity range,
  virtual void PrintInfo() override;

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*) override;

private:
  
  // hide assignment operator as private 
  G4PolarizedGammaConversion& 
    operator=(const G4PolarizedGammaConversion &right) = delete;
  G4PolarizedGammaConversion(const G4PolarizedGammaConversion& ) = delete;
     
  G4bool          isInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4PolarizedGammaConversion::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Gamma::Gamma());
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
#endif
 
