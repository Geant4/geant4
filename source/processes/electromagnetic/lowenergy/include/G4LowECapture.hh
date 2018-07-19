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
// $Id: G4LowECapture.hh,v 1.1 2010-08-31 11:23:58 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4LowECapture
//
// Description: The process to kill low-energy particles
//
// Author:      V.Ivanchenko 12 May 2015
//
//----------------------------------------------------------------------------
//
// Class description:
//
// G4LowECapture allows to remove unwanted particles from simulation in 
// order to improve CPU performance. There are two parameters:
//                 
// 1) low energy threshold for kinetic energy (default 0)
// 2) the vector of names of G4Region where process is active
// 
//
// If a track is killed then energy deposition is added to the step 
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef LowECapture_h
#define LowECapture_h 1

#include "G4VDiscreteProcess.hh"
#include "globals.hh"

class G4Region;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4LowECapture : public G4VDiscreteProcess
{
public:

  G4LowECapture(G4double ekinlimit = 0.0);

  virtual ~G4LowECapture();

  void SetKinEnergyLimit(G4double);

  void AddRegion(const G4String&);

  virtual void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
							G4double,
							G4ForceCondition*);

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:

  virtual G4double GetMeanFreePath(const G4Track&, G4double,G4ForceCondition*);

private:

  // hide assignment operator as private
  G4LowECapture(const G4LowECapture&);
  G4LowECapture& operator = (const G4LowECapture &right);

  G4double kinEnergyThreshold;
  G4bool   isIon;
  G4int    nRegions;
  std::vector<G4String> regionName;
  std::vector<const G4Region*> region;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
