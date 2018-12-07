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
// -------------------------------------------------------------------
//
// Author: A. Forti
//
// History:
// -----------
// 02 Mar 1999  A. Forti        1st implementation
// 27 Sep 2001  V. Ivanchenko   Major revision according to a design iteration
// 10 Oct 2001  M.G. Pia        Revision to improve code quality and 
//                              consistency with design
// 18 Oct 2001  M.G. Pia        Revision to improve code quality and 
//                              consistency with design
// 29 Nov 2001  V.Ivanchenko    New parametrisation of EEDL data
// 31 May 2002  V.Ivanchenko    Add Auger flag
// 06 Feb 2003  V.Ivanchenko    Change signature of deexcitation for cut per region
//
// -------------------------------------------------------------------

// Class description:
// Low Energy electromagnetic process, electron Ionisation
// based on the data of the EEDL database. Details are described
// in the Physics Reference Manual.
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// --------------------------------------------------------------

#ifndef G4RDLOWENERGYIONISATION_HH
#define G4RDLOWENERGYIONISATION_HH 1

#include "G4eLowEnergyLoss.hh"
#include "G4RDAtomicDeexcitation.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4RDVDataSetAlgorithm;
class G4ParticleChange;
class G4RDVEnergySpectrum;
class G4RDVCrossSectionHandler;
class G4RDShellVacancy;
class G4RDVEMDataSet;

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

  void ActivateAuger(G4bool val);
    
protected:
 
  G4double GetMeanFreePath(const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition );

protected:

  virtual std::vector<G4DynamicParticle*>* DeexciteAtom(const G4MaterialCutsCouple* couple,
							  G4double incidentEnergy,
							  G4double eLoss);

private:

  // Hide copy constructor and assignment operator as private 
  G4LowEnergyIonisation(const G4LowEnergyIonisation& );
  G4LowEnergyIonisation& operator = (const G4LowEnergyIonisation& right);
  
  void BuildLossTable(const G4ParticleDefinition& ParticleType);
  
  G4RDVCrossSectionHandler* crossSectionHandler;
  G4RDVEMDataSet* theMeanFreePath;
  G4RDVEnergySpectrum* energySpectrum;

  // Lower limit for generation of gamma in this model
  G4DataVector cutForDelta;
  G4double cutForPhotons;
  G4double cutForElectrons;
  G4RDAtomicDeexcitation deexcitationManager;
  G4RDShellVacancy* shellVacancy;
  
};

#endif
 










