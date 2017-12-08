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
// $Id: G4MuMultipleScattering.hh 107056 2017-11-01 14:52:32Z gcosmo $
//
// -----------------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4MuMultipleScattering
//
// Author:        V.Ivanchenko using Laszlo Urban original code
//
// Creation date: 24.10.2007 cloned from G4MultipleScattering by VI
// 
// Modifications:
// 20.03.07 Remove local parameter skin (V.Ivanchenko) 
//
//
//------------------------------------------------------------------------------

// class description
//
//  The class simulates the multiple scattering for any kind
//  of charged particle.
//
// class description - end

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4MuMultipleScattering_h
#define G4MuMultipleScattering_h 1

#include "G4VMultipleScattering.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MuMultipleScattering : public G4VMultipleScattering
{
public:    // with description

  explicit G4MuMultipleScattering(const G4String& processName="muMsc");

  virtual ~G4MuMultipleScattering();

  // returns true for charged particles, false otherwise
  G4bool IsApplicable (const G4ParticleDefinition& p) override;

  // print description in html
  virtual void ProcessDescription(std::ostream&) const override;

protected:

  // Print out of the class parameters
  virtual void StreamProcessInfo(std::ostream& outFile,
                             G4String endOfLine=G4String("\n")) const override;

  // This function initialise models
  void InitialiseProcess(const G4ParticleDefinition*) override;

private:        // data members

  //  hide assignment operator
  G4MuMultipleScattering & 
    operator=(const  G4MuMultipleScattering &right) = delete;
  G4MuMultipleScattering(const  G4MuMultipleScattering&) = delete;

  G4bool   isInitialized;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
