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
// $Id: G4PenelopeRayleigh.hh,v 1.1 2002-12-06 16:25:21 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 01 Dec 2002   L.Pandola  1st implementation
// -------------------------------------------------------------------
// Class description:
// Low Energy Electromagnetic process, Rayleigh effect
// Penelope model
// -------------------------------------------------------------------

#ifndef G4PENELOPERAYLEIGH_HH
#define G4PENELOPERAYLEIGH_HH 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4PenelopeIntegrator.hh"
class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VEMDataSet;
class G4Material;
class G4DataVector;

class G4PenelopeRayleigh : public G4VDiscreteProcess {

public:
  
  G4PenelopeRayleigh(const G4String& processName ="PenRayleigh");
  
  ~G4PenelopeRayleigh();

  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void BuildPhysicsTable(const G4ParticleDefinition& photon);
  
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);

  G4double MolecularFormFactor(G4double x); 
  G4double DifferentialCrossSection (G4double e);
 
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
  G4PenelopeRayleigh& operator=(const G4PenelopeRayleigh &right);
  G4PenelopeRayleigh(const G4PenelopeRayleigh& );
  
  G4double lowEnergyLimit;  // low energy limit  applied to the process
  G4double highEnergyLimit; // high energy limit applied to the process

  void InizialiseSampling();

  G4Material* material;
  G4DataVector* samplingFunction_x;
  G4DataVector* samplingFunction_y;

  G4double samplingConstant;
  G4double facte; //cross section factor

  G4VEMDataSet* meanFreePathTable;
  
  const G4int nBins;
  const G4double intrinsicLowEnergyLimit; // intrinsic validity range
  const G4double intrinsicHighEnergyLimit;
};

#endif

 












