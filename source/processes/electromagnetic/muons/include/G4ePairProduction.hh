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
// $Id: G4ePairProduction.hh 72895 2013-08-09 11:34:16Z vnivanch $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4ePairProduction
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 17.03.2016
//
// Modifications:
//
// Class Description:
//
// This class manages the PairProduction process for e+-.
// it inherites from G4VEnergyLossProcess.
//

// -------------------------------------------------------------------
//

#ifndef G4ePairProduction_h
#define G4ePairProduction_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4VEnergyLossProcess.hh"

class G4ePairProduction : public G4VEnergyLossProcess
{
public:

  explicit G4ePairProduction(const G4String& processName = "ePairProd");

  virtual ~G4ePairProduction();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p) override;

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
				    const G4Material*, G4double cut) override;

  inline void SetLowestKineticEnergy(G4double e);

  // print description in html
  virtual void ProcessDescription(std::ostream&) const override;

protected:

  // Print out of the class parameters
  virtual void StreamProcessInfo(std::ostream& outFile,
                             G4String endOfLine=G4String("\n")) const override;

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
					   const G4ParticleDefinition*) override;

private:

  G4ePairProduction & operator=(const G4ePairProduction &right) = delete;
  G4ePairProduction(const G4ePairProduction&) = delete;

protected:

  const G4ParticleDefinition* theParticle;
  G4double                    lowestKinEnergy;
  G4bool                      isInitialised;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4ePairProduction::SetLowestKineticEnergy(G4double e) 
{
  lowestKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
