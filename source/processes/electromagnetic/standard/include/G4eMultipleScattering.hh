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
// $Id: G4eMultipleScattering.hh,v 1.4 2010-10-26 10:39:02 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4eMultipleScattering
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

#ifndef G4eMultipleScattering_h
#define G4eMultipleScattering_h 1

#include "G4VMultipleScattering.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4eMultipleScattering : public G4VMultipleScattering

{
public:    // with description

  G4eMultipleScattering(const G4String& processName = "msc");

  virtual ~G4eMultipleScattering();

  // returns true for charged particles, false otherwise
  G4bool IsApplicable (const G4ParticleDefinition& p);

  // Print few lines of informations about the process: validity range,
  void PrintInfo();

  // geom. step length distribution should be sampled or not
  //void Setsamplez(G4bool value) { samplez = value;};

  // to reduce the energy/step dependence
  //void Setdtrl(G4double value) { dtrl = value;};

  // 'soften' step limitation above lambdalimit
  //void SetLambdalimit(G4double value) { lambdalimit = value;};

protected:

  // This function initialise models
  void InitialiseProcess(const G4ParticleDefinition*);

private:        // data members

  //  G4double lambdalimit;
  //G4double dtrl;

  //G4bool   samplez;
  G4bool   isInitialized;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
