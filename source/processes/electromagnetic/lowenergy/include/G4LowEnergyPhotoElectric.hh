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
// $Id: G4LowEnergyPhotoElectric.hh,v 1.22 2001-09-21 08:43:32 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: A. Forti
//         Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 02 Mar  1999   A. Forti   1st implementation
// 12 Aug  2001   MGP        Major revision according to a design iteration
// 16 Sept 2001   E. Guardincerri Added fluorescence generation
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
#include "G4ShellData.hh"
#include "G4AtomicDeexcitation.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VEMDataSet;
class G4CrossSectionHandler;
class G4VDataSetAlgorithm;

class G4LowEnergyPhotoElectric : public G4VDiscreteProcess {
  
public:
  
  G4LowEnergyPhotoElectric(const G4String& processName ="LowEnPhotoElec");
  
  ~G4LowEnergyPhotoElectric();

  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void BuildPhysicsTable(const G4ParticleDefinition& photon);
 
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);

  void SetCutForLowEnSecPhotons(G4double);

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

  G4VDataSetAlgorithm* scatterInterpolation;

  G4VEMDataSet* meanFreePathTable;

  G4CrossSectionHandler* crossSectionHandler;
  G4CrossSectionHandler* shellCrossSectionHandler;

  const G4double intrinsicLowEnergyLimit; // intrinsic validity range
  const G4double intrinsicHighEnergyLimit;

  G4double cutForLowEnergySecondaryPhotons;

  G4ShellData shellData;
  G4AtomicDeexcitation deexcitationManager;
};

#endif



