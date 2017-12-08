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
// $Id: G4CoulombScattering.hh 106717 2017-10-20 09:41:27Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4CoulombScattering
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 12.03.2006
//
// Modifications:
//
// Class Description:
//
// This class manages the process of Coulomb elastic scattering
//

// -------------------------------------------------------------------
//

#ifndef G4CoulombScattering_h
#define G4CoulombScattering_h 1

#include "G4VEmProcess.hh"
#include "G4VEmModel.hh"

class G4CoulombScattering : public G4VEmProcess
{

public:

  explicit G4CoulombScattering(const G4String& name = "CoulombScat");

  virtual ~G4CoulombScattering();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) final;

  // print documentation in html format
  virtual void ProcessDescription(std::ostream&) const override;

protected:

  // Print out of the class parameters
  virtual void StreamProcessInfo(std::ostream& outFile,
                             G4String endOfLine=G4String("\n")) const override;

  virtual void InitialiseProcess(const G4ParticleDefinition*) override;

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
                                    const G4Material*) final;

private:

 // hide assignment operator
  G4CoulombScattering & operator=(const G4CoulombScattering &right) = delete;
  G4CoulombScattering(const G4CoulombScattering&) = delete;
  
  G4double q2Max;
  G4bool isInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
