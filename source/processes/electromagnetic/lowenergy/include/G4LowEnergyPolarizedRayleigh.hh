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
// $Id: G4LowEnergyPolarizedRayleigh.hh,v 1.1 2003-05-16 08:58:39 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
// Author: A. Forti
//         Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//         Gerardo O. Depaola 
//
// History:
// -----------
//
// 25 April 2003   G.Depaola & F.Longo 1st implementation
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Rayleigh effect & Polarization

// -------------------------------------------------------------------

#ifndef G4LOWENERGYPOLARIZEDRAYLEIGH_HH
#define G4LOWENERGYPOLARIZEDRAYLEIGH_HH 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VEMDataSet;
class G4VCrossSectionHandler;
class G4VRangeTest;

class G4LowEnergyPolarizedRayleigh : public G4VDiscreteProcess {

public:
  
  G4LowEnergyPolarizedRayleigh(const G4String& processName ="polarLowEnRayleigh");
  
  ~G4LowEnergyPolarizedRayleigh();

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
  G4LowEnergyPolarizedRayleigh& operator=(const G4LowEnergyPolarizedRayleigh &right);
  G4LowEnergyPolarizedRayleigh(const G4LowEnergyPolarizedRayleigh& );
  
  G4double lowEnergyLimit;  // low energy limit  applied to the process
  G4double highEnergyLimit; // high energy limit applied to the process

  G4VEMDataSet* meanFreePathTable;
  G4VEMDataSet* formFactorData;

  G4VCrossSectionHandler* crossSectionHandler;
  G4VRangeTest* rangeTest;

  const G4double intrinsicLowEnergyLimit; // intrinsic validity range
  const G4double intrinsicHighEnergyLimit;
  G4ThreeVector GetRandomPolarization(G4ThreeVector& direction0); // Random Polarization
  G4ThreeVector GetPerpendicularPolarization(const G4ThreeVector& direction0, const G4ThreeVector& polarization0) const;
  
  G4ThreeVector SetPerpendicularVector(G4ThreeVector& a); 
  G4ThreeVector SetNewPolarization(G4double sinSqrTheta, 
				   G4double phi, G4double cosTheta);
  G4double SetPhi(G4double);
  
  void SystemOfRefChange(G4ThreeVector& direction0, G4ThreeVector& direction1, 
			 G4ThreeVector& polarization0, G4ThreeVector& polarization1);

};

#endif

 












