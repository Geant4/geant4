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
// $Id: G4LowEnergyOldCompton.hh,v 1.1 2001-09-12 17:01:14 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file --- Copyright CERN 1995
//      CERN Geneva Switzerland
//
//      ------------ G4LowEnergyCompton physics process ------
//                   by A.Forti 1999/03/02
//
// Class description:
// Low Energy electromagnetic process, Compton
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// ************************************************************

#ifndef G4LowEnergyCompton_h
#define G4LowEnergyCompton_h 

// Base Class Headers
#include "G4VDiscreteProcess.hh"

// Contained Variables Headers
#include "G4LowEnergyUtilities.hh"
#include "G4Gamma.hh"

class G4LowEnergyCompton : public G4VDiscreteProcess{

private: 

  // hide assignment operator as private 
  G4LowEnergyCompton& operator=(const G4LowEnergyCompton &right);
  G4LowEnergyCompton(const G4LowEnergyCompton& );
 
public:
  
  G4LowEnergyCompton(const G4String& processName ="LowEnCompton");
  
  ~G4LowEnergyCompton();

  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void BuildPhysicsTable(const G4ParticleDefinition& GammaType);
 
  G4double GetMeanFreePath(const G4Track& aTrack, 
			   G4double previousStepSize, 
			   G4ForceCondition* condition);

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
  
protected:

  void BuildScatteringFunctionTable();
  void BuildCrossSectionTable();
  void BuildMeanFreePathTable();
  void BuildZVec();

private:

  G4Element* SelectRandomAtom(const G4DynamicParticle*, G4Material*);
  
  G4SecondLevel* theCrossSectionTable;
  G4SecondLevel* theScatteringFunctionTable;
  G4PhysicsTable* theMeanFreePathTable;
  G4DataVector* ZNumVec;

  G4double lowestEnergyLimit; // low  energy limit of the crosssection data 
  G4double highestEnergyLimit; // high energy limit of the crosssection data
  G4int numbBinTable; // number of bins in the data  tables

  G4LowEnergyUtilities util;

  G4double meanFreePath; // actual Mean Free Path (current medium)
};

#include "G4LowEnergyCompton.icc"

#endif

