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
// -------------------------------------------------------------------
// $Id: G4LowEnergyIonisation.hh,v 1.33 2001-11-29 19:01:45 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: A. Forti
//
// History:
// -----------
// 02 Mar 1999  A. Forti        1st implementation
// 27 Sep 2001  V. Ivanchenko   Major revision according to a design iteration
// 10 Oct 2001  M.G. Pia        Revision to improve code quality and consistency with design
// 18 Oct 2001  M.G. Pia        Revision to improve code quality and consistency with design
//
// -------------------------------------------------------------------

// Class description:
// Low Energy electromagnetic process, electron Ionisation
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// --------------------------------------------------------------

#ifndef G4lOWENERGYIONISATION_HH
#define G4LOWENERGYIONISATION_HH 1

#include "G4eLowEnergyLoss.hh"
#include "G4AtomicDeexcitation.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VDataSetAlgorithm;
class G4ParticleChange;
class G4VEnergySpectrum;
class G4VCrossSectionHandler;
class G4ShellVacancy;
class G4VEMDataSet;

class G4LowEnergyIonisation : public G4eLowEnergyLoss
{ 
public:
 
  G4LowEnergyIonisation(const G4String& processName = "LowEnergyIoni");
  
  ~G4LowEnergyIonisation();
  
  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void PrintInfoDefinition();
  
  void BuildPhysicsTable(const G4ParticleDefinition& ParticleType);
  
  G4VParticleChange* PostStepDoIt(const G4Track& track,         
				  const G4Step& step);                 
 
  void SetCutForLowEnSecPhotons(G4double cut);

  void SetCutForLowEnSecElectrons(G4double cut);
    
 
  G4double GetMeanFreePath(const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition );
protected:

  virtual G4std::vector<G4DynamicParticle*>* DeexciteAtom(const G4Material* material,
							  G4double incidentEnergy,
							  G4double eLoss);

private:

  // Hide copy constructor and assignment operator as private 
  G4LowEnergyIonisation(const G4LowEnergyIonisation& );
  G4LowEnergyIonisation& operator = (const G4LowEnergyIonisation& right);
  
  void BuildLossTable(const G4ParticleDefinition& ParticleType);
  
  G4VCrossSectionHandler* crossSectionHandler;
  G4VEMDataSet* theMeanFreePath;
  G4VEnergySpectrum* energySpectrum;

  // Lower limit for generation of gamma in this model
  G4DataVector cutForDelta;
  G4double cutForPhotons;
  G4double cutForElectrons;
  G4AtomicDeexcitation deexcitationManager;
  G4ShellVacancy* shellVacancy;
  
};

#endif
 










