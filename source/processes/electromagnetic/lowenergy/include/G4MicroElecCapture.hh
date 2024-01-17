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
// $Id: G4ElectronCapture.hh,v 1.1 2010-08-31 11:23:58 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4MicroElecCapture
//
// Description: The process to kill e- to save CPU
//
// Author:      V.Ivanchenko 31 August 2010 modified and adapted to MicorElec by C. Inguimbert 321/01/2022
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

#ifndef MicroElecCapture_h
#define MicroElecCapture_h 1

#include "G4VDiscreteProcess.hh"
#include "G4MicroElecMaterialStructure.hh"
#include "globals.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleChangeForGamma.hh"

class G4Region;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MicroElecCapture : public G4VDiscreteProcess
{
public:

  G4MicroElecCapture(const G4String& regName, G4double ekinlimit);

  virtual ~G4MicroElecCapture();

  void SetKinEnergyLimit(G4double);

  void BuildPhysicsTable(const G4ParticleDefinition&) override;

  G4bool IsApplicable(const G4ParticleDefinition&) override;

  void Initialise();
  
  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

  G4MicroElecCapture(const G4MicroElecCapture&) = delete;
  G4MicroElecCapture& operator = (const G4MicroElecCapture &right) = delete;


protected:

  G4double GetMeanFreePath(const G4Track&, G4double,G4ForceCondition*) override;

private:

  G4double G_Lindhard_Rob(G4double , G4int , G4int , G4int , G4int );
  
  typedef std::map<G4String, G4MicroElecMaterialStructure*, std::less<G4String> > WorkFunctionTable;
  WorkFunctionTable tableWF; //Table of all materials simulated

  G4bool isInitialised;
  G4double kinEnergyThreshold;
  G4String regionName;
  G4Region* region;
  G4ParticleChangeForGamma fParticleChange;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

