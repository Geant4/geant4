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
// -------------------------------------------------------------------
// $Id: G4PenelopeBremsstrahlung.hh,v 1.7 2003/11/11 17:18:30 pandola Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Author: L.Pandola
//
// History:
// -----------
// 03 Feb 2003  L. Pandola       1st implementation
// 18 Mar 2003  L. Pandola       positrons added
// 23 May 2003  L. Pandola       ebins added as private member
// 25 May 2003  MGP              ebins removed from data members
// 11 Nov 2003  L. Pandola       Code review: use std::map for angular data
// Class description:
// Penelope electromagnetic process, electron and positron Bremsstrahlung
// --------------------------------------------------------------


#ifndef G4PENELOPEBREMSSTRAHLUNG_HH
#define G4PENELOPEBREMSSTRAHLUNG_HH 1

#include "G4eLowEnergyLoss.hh"
#include "G4DataVector.hh"
#include "globals.hh"
#include "G4PenelopeBremsstrahlungAngular.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VEMDataSet;
class G4VEnergySpectrum;
class G4VCrossSectionHandler;
class G4PenelopeBremsstrahlung : public G4eLowEnergyLoss
{ 
public:
 
  G4PenelopeBremsstrahlung(const G4String& processName = "PenelopeBrem");
  
  ~G4PenelopeBremsstrahlung();
  
  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void BuildPhysicsTable(const G4ParticleDefinition& particleType);
  
  G4VParticleChange* PostStepDoIt(const G4Track& track,         
				  const G4Step& step);                 
 
  void SetCutForLowEnSecPhotons(G4double cut);
  
  void PrintInfoDefinition();

  //For testing purpose only
  G4double DumpMeanFreePath(const G4Track& aTrack,
                            G4double previousStepSize,
                            G4ForceCondition* condition)
  { return GetMeanFreePath(aTrack, previousStepSize, condition); }
  
protected:

  G4double GetMeanFreePath(const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition );
 
private:

  // Hide copy constructor and assignment operator as private 
  G4PenelopeBremsstrahlung(const G4PenelopeBremsstrahlung& );
  G4PenelopeBremsstrahlung& operator = (const G4PenelopeBremsstrahlung& right);
  
  void BuildLossTable(const G4ParticleDefinition& ParticleType);
  
  void LoadAngularData();

  G4VCrossSectionHandler* crossSectionHandler;
  G4VEMDataSet* theMeanFreePath;
  G4VEnergySpectrum* energySpectrum;

  // Map to the objects containing tha angular data 
  std::map<G4int,G4PenelopeBremsstrahlungAngular*> *angularData;

 
  // Lower limit for generation of gamma in this model
  G4DataVector cutForSecondaryPhotons;
  G4double cutForPhotons;

};


  
#endif
 










