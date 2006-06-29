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
// $Id: G4PenelopePhotoElectric.hh,v 1.4 2006-06-29 19:36:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: L.Pandola
//
// History:
// -----------
// 07 Jan  2003   L.Pandola   1st implementation
// 31 May  2005   L.Pandola   Added method to sample angular distribution
//                            of the final electron
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic process, Photoelectric effect with
// Penelope model and tables
// -------------------------------------------------------------------

#ifndef G4PENELOPEPHOTOELECTRIC_HH
#define G4PENELOPEPHOTOELECTRIC_HH 1

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

class G4PenelopePhotoElectric : public G4VDiscreteProcess {
  
public:
  
  G4PenelopePhotoElectric(const G4String& processName ="PenPhotoElec");
  
  ~G4PenelopePhotoElectric();

  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void BuildPhysicsTable(const G4ParticleDefinition& photon);
 
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);

  void SetCutForLowEnSecPhotons(G4double);

  void SetCutForLowEnSecElectrons(G4double);

  void ActivateAuger(G4bool);

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
  G4PenelopePhotoElectric& operator=(const G4PenelopePhotoElectric& right);
  G4PenelopePhotoElectric(const G4PenelopePhotoElectric& );

  G4double SampleElectronDirection(G4double);

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
};

#endif



