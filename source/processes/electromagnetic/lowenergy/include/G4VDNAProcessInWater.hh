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
// -------------------------------------------------------------------
// $Id: G4VDNAProcessInWater.hh,v 1.7 2007-10-08 09:18:43 sincerti Exp $
// -------------------------------------------------------------------
//

#ifndef G4VDNAProcessInWater_HH
#define G4VDNAProcessInWater_HH 1
 
#include "G4VLowEnergyTestableDiscreteProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

template<typename TotalCrossSectionPolicy, typename FinalStatesPolicy> 
class G4VDNAProcessInWater:
public G4VLowEnergyTestableDiscreteProcess, 
public TotalCrossSectionPolicy, 
public FinalStatesPolicy
{
 public:

  G4VDNAProcessInWater(const G4String& name) : G4VLowEnergyTestableDiscreteProcess(name) {}

  virtual ~G4VDNAProcessInWater() {}
 
  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);

  virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleDefinition);

 protected:

  virtual G4double GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize, G4ForceCondition* condition);

 private:

  // Hides default constructor and assignment operator as private 
  G4VDNAProcessInWater();
  G4VDNAProcessInWater& operator=(const G4VDNAProcessInWater & right);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VDNAProcessInWater.icc"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif 
