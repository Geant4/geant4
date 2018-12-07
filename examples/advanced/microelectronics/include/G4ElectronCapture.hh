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
//---------------------------------------------------------------------------
//
// ClassName:   G4ElectronCapture
//
// Description: The process to kill e- to save CPU
//
// Author:      V.Ivanchenko 31 August 2010
//
//----------------------------------------------------------------------------
//
// Class description:
//
// G4ElectronCapture allows to remove unwanted e- from simulation in 
// order to improve CPU performance. There are two parameters:
//                 
// 1) low energy threshold for e- kinetic energy (default 0)
// 2) the name of G4Region where process is active
// 
//
// If an electron track is killed then energy deposition is added to the step 
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ElectronCapture_h
#define ElectronCapture_h 1

#include "G4VDiscreteProcess.hh"
#include "globals.hh"
#include "G4ParticleChangeForGamma.hh"

class G4Region;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ElectronCapture : public G4VDiscreteProcess
{
public:

  G4ElectronCapture(const G4String& regName, G4double ekinlimit);

  virtual ~G4ElectronCapture();

  void SetKinEnergyLimit(G4double);

  virtual void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  virtual G4double PostStepGetPhysicalInteractionLength( const G4Track& track,
							 G4double previousStepSize,
							 G4ForceCondition* condition);

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:

  virtual G4double GetMeanFreePath(const G4Track&, G4double,G4ForceCondition*);

private:

  // hide assignment operator as private
  G4ElectronCapture(const G4ElectronCapture&);
  G4ElectronCapture& operator = (const G4ElectronCapture &right);

  G4double kinEnergyThreshold;
  G4String regionName;
  G4Region* region;
  G4ParticleChangeForGamma fParticleChange;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

