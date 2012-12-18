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
// $Id: G4hBremsstrahlung.hh,v 1.2 2009-02-20 16:38:33 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4hBremsstrahlung
//
// Author:        Vladimir Ivanchenko on base of model for muons
//
// Creation date: 01.03.2008
//
// Modifications:
//
//
// Class Description:
//
// This class manages the Bremsstrahlung process for hadrons
// it inherites from G4VContinuousDiscreteProcess via G4VEnergyLossProcess.
//

// -------------------------------------------------------------------
//

#ifndef G4hBremsstrahlung_h
#define G4hBremsstrahlung_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4VEmModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleDefinition;

class G4hBremsstrahlung : public G4VEnergyLossProcess

{
public:

  G4hBremsstrahlung(const G4String& processName = "hBrems");

  virtual ~G4hBremsstrahlung();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p);

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
				    const G4Material*, 
				    G4double cut);

  // Print out of the class parameters
  virtual void PrintInfo();

protected:

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
					   const G4ParticleDefinition*);

private:

  G4hBremsstrahlung & operator=(const G4hBremsstrahlung &right);
  G4hBremsstrahlung(const G4hBremsstrahlung&);

  const G4ParticleDefinition* theParticle;
  const G4ParticleDefinition* theBaseParticle;

  G4double  lowestKinEnergy;
  G4bool    isInitialised;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
