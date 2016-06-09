//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4LowEnergyPhotoElectric.hh,v 1.26 2004/05/12 09:20:03 silvarod Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

#ifndef G4LOWENERGYPHOTOELECTRIC_HH
#define G4LOWENERGYPHOTOELECTRIC_HH 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4AtomicDeexcitation.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VEMDataSet;
class G4VCrossSectionHandler;
class G4VRangeTest;
class G4VPhotoElectricAngularDistribution;

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

  void SetAngularGenerator(G4VPhotoElectricAngularDistribution* distribution);

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

  G4VEMDataSet* meanFreePathTable;

  G4VCrossSectionHandler* crossSectionHandler;
  G4VCrossSectionHandler* shellCrossSectionHandler;

  G4VRangeTest* rangeTest;

  const G4double intrinsicLowEnergyLimit; // intrinsic validity range
  const G4double intrinsicHighEnergyLimit;

  G4double cutForLowEnergySecondaryPhotons;
  G4double cutForLowEnergySecondaryElectrons;

  G4AtomicDeexcitation deexcitationManager;

  G4VPhotoElectricAngularDistribution* ElectronAngularGenerator;
  G4String generatorName;
};

#endif



