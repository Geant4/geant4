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
// $Id: G4eAdjointMultipleScattering.hh 96934 2016-05-18 09:10:41Z gcosmo $
//
// -----------------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4eAdjointMultipleScattering
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 10 March 2001
// 
// Modifications:
//
//
//------------------------------------------------------------------------------
//

// class description
//
//  The class simulates the multiple scattering for e+ and e-
//
// class description - end

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4eAdjointMultipleScattering_h
#define G4eAdjointMultipleScattering_h 1

#include "G4VMultipleScattering.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4eAdjointMultipleScattering : public G4VMultipleScattering

{
public:    // with description

  explicit G4eAdjointMultipleScattering(const G4String& processName = "msc");

  virtual ~G4eAdjointMultipleScattering();

  // This is called in the beginning of tracking for a new track
  void StartTracking(G4Track*) override;

  // returns true for charged particles, false otherwise
  G4bool IsApplicable (const G4ParticleDefinition& p) final;

  // Print few lines of informations about the process: validity range,
  void PrintInfo() override;

protected:

  // This function initialise models
  void InitialiseProcess(const G4ParticleDefinition*) override;

private:        // data members

  G4bool   isInitialized;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
