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
// $Id: G4PenelopeBremsstrahlung.hh,v 1.2 2003/03/19 10:27:44 pandola Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// Author: L.Pandola
//
// History:
// -----------
// 03 Feb 2003  L. Pandola       1st implementation
// 18 Mar 2003  L. Pandola       positrons added
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
  typedef G4std::vector<G4PenelopeBremsstrahlungAngular*> G4AngularData;
  //vector of pointers to the angular factors of the elements in each material

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
  G4std::vector<G4AngularData*> materialAngularData;
  //Vector of pointers to the vectors containing tha angular data 
  //one element=one material
 
  // Lower limit for generation of gamma in this model
  G4DataVector cutForSecondaryPhotons;
  G4double cutForPhotons;
};


  
#endif
 










