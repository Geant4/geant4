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
// Author: A. Forti
//         Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 02 Mar 1999   A. Forti   1st implementation
//  1 Aug 2001   MGP        Major revision according to a design iteration
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics, Compton Scattering
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4RDLOWENERGYCOMPTON_HH
#define G4RDLOWENERGYCOMPTON_HH 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4RDShellData.hh"
#include "G4RDDopplerProfile.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4RDVEMDataSet;
class G4RDVCrossSectionHandler;
class G4RDVRangeTest;

class G4LowEnergyCompton : public G4VDiscreteProcess {

public:
  
  G4LowEnergyCompton(const G4String& processName ="LowEnCompton");
  
  ~G4LowEnergyCompton();

  G4bool IsApplicable(const G4ParticleDefinition& definition);
  
  void BuildPhysicsTable(const G4ParticleDefinition& definition);
 
  G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);
 
  // For testing purpose only
  G4double DumpMeanFreePath(const G4Track& track, 
			    G4double previousStepSize, 
			    G4ForceCondition* condition) 
  { return GetMeanFreePath(track, previousStepSize, condition); }

protected:

  G4double GetMeanFreePath(const G4Track& track, 
			   G4double previousStepSize, 
			   G4ForceCondition* condition);

private: 

  // Hide copy constructor and assignment operator as private 
  G4LowEnergyCompton& operator=(const G4LowEnergyCompton& right);
  G4LowEnergyCompton(const G4LowEnergyCompton& );

  G4double lowEnergyLimit;  // low energy limit  applied to the process
  G4double highEnergyLimit; // high energy limit applied to the process

  G4RDVEMDataSet* meanFreePathTable;
  G4RDVEMDataSet* scatterFunctionData;

  G4RDVCrossSectionHandler* crossSectionHandler;

  G4RDVRangeTest* rangeTest;

  const G4double intrinsicLowEnergyLimit; // intrinsic validity range
  const G4double intrinsicHighEnergyLimit;

  G4RDShellData shellData;
  G4RDDopplerProfile profileData;
};

#endif

