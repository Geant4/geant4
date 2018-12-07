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
//
// Author: A. Forti
//         Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 02 Mar  1999   A. Forti   1st implementation
// 12 Aug  2001   MGP        Major revision according to a design iteration
// 16 Sept 2001   E. Guardincerri Added fluorescence generation
// 31 May  2002   V.Ivanchenko Add cut on Auger electrons
// 10 May  2004   P. Rodrigues Changes to accommodate new angular generators
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic process, Photoelectric effect
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4RDLOWENERGYPHOTOELECTRIC_HH
#define G4RDLOWENERGYPHOTOELECTRIC_HH 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4RDAtomicDeexcitation.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4RDVEMDataSet;
class G4RDVCrossSectionHandler;
class G4RDVRangeTest;
class G4RDVPhotoElectricAngularDistribution;

class G4LowEnergyPhotoElectric : public G4VDiscreteProcess {
  
public:
  
  G4LowEnergyPhotoElectric(const G4String& processName ="LowEnPhotoElec");
  
  ~G4LowEnergyPhotoElectric();

  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void BuildPhysicsTable(const G4ParticleDefinition& photon);
 
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);

  void SetCutForLowEnSecPhotons(G4double);

  void SetCutForLowEnSecElectrons(G4double);

  void ActivateAuger(G4bool);

  void SetAngularGenerator(G4RDVPhotoElectricAngularDistribution* distribution);

  void SetAngularGenerator(const G4String& name);


  // For testing purpose only
  G4double DumpMeanFreePath(const G4Track& aTrack, 
			    G4double previousStepSize, 
			    G4ForceCondition* condition) 
  { return GetMeanFreePath(aTrack, previousStepSize, condition); }

protected:  

  G4double GetMeanFreePath(const G4Track& aTrack, 
			   G4double previousStepSize, 
			   G4ForceCondition* condition);

private:

 // Hide copy constructor and assignment operator as private 
  G4LowEnergyPhotoElectric& operator=(const G4LowEnergyPhotoElectric& right);
  G4LowEnergyPhotoElectric(const G4LowEnergyPhotoElectric& );

  G4double lowEnergyLimit;  // low energy limit  applied to the process
  G4double highEnergyLimit; // high energy limit applied to the process

  G4RDVEMDataSet* meanFreePathTable;

  G4RDVCrossSectionHandler* crossSectionHandler;
  G4RDVCrossSectionHandler* shellCrossSectionHandler;

  G4RDVRangeTest* rangeTest;

  const G4double intrinsicLowEnergyLimit; // intrinsic validity range
  const G4double intrinsicHighEnergyLimit;

  G4double cutForLowEnergySecondaryPhotons;
  G4double cutForLowEnergySecondaryElectrons;

  G4RDAtomicDeexcitation deexcitationManager;

  G4RDVPhotoElectricAngularDistribution* ElectronAngularGenerator;
  G4String generatorName;
};

#endif



