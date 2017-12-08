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
// $Id: G4hBremsstrahlung.hh 106715 2017-10-20 09:39:06Z gcosmo $
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
#include "G4MuBremsstrahlung.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4hBremsstrahlung : public G4MuBremsstrahlung

{
public:

  explicit G4hBremsstrahlung(const G4String& processName = "hBrems");

  virtual ~G4hBremsstrahlung();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) override;

  // print description in html
  virtual void ProcessDescription(std::ostream&) const override;

protected:

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
					   const G4ParticleDefinition*) override;

private:

  G4hBremsstrahlung & operator=(const G4hBremsstrahlung &right) = delete;
  G4hBremsstrahlung(const G4hBremsstrahlung&) = delete;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
