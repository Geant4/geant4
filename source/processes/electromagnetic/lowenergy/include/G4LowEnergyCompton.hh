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

// $Id: G4LowEnergyCompton.hh,v 1.17 2001-09-10 18:05:16 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#ifndef G4LOWENERGYCOMPTON_HH
#define G4LOWENERGYCOMPTON_HH 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VEMDataSet;
class G4CrossSectionHandler;
class G4VDataSetAlgorithm;

class G4LowEnergyCompton : public G4VDiscreteProcess {

public:
  
  G4LowEnergyCompton(const G4String& processName ="LowEnCompton");
  
  ~G4LowEnergyCompton();

  G4bool IsApplicable(const G4ParticleDefinition& definition);
  
  void BuildPhysicsTable(const G4ParticleDefinition& photon);
 
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
 
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
  G4LowEnergyCompton& operator=(const G4LowEnergyCompton& right);
  G4LowEnergyCompton(const G4LowEnergyCompton& );

  G4double lowEnergyLimit;  // low energy limit  applied to the process
  G4double highEnergyLimit; // high energy limit applied to the process

  G4VDataSetAlgorithm* scatterInterpolation;

  G4VEMDataSet* meanFreePathTable;
  G4VEMDataSet* scatterFunctionData;

  G4CrossSectionHandler* crossSectionHandler;

  const G4double intrinsicLowEnergyLimit; // intrinsic validity range
  const G4double intrinsicHighEnergyLimit;

};

#endif

