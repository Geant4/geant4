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
// $Id: G4LowEnergyPolarizedCompton.hh,v 1.6 2001-09-10 18:05:16 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//

// --------- G4LowEnergyPolarizedCompton class -----
//
//           by G.Depaola & F.Longo (21 may 2001)
// 24 May 2001 - MGP      Modified to inherit from G4VDiscreteProcess
// 25 May 2001 - MGP      Added protections to avoid crashes
//
// Class description:
// Low Energy electromagnetic process, Polarised Compton scattering
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// ------------------------------------------------------------

#ifndef G4LowEnergyPolarizedCompton_h
#define G4LowEnergyPolarizedCompton_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4LowEnergyUtilities.hh"
//#include "G4VParticleChange.hh"

class G4SecondLevel;
class G4PhysicsTable;
class G4DataVector;
class G4ParticleDefinition;
class G4VParticleChange;


class G4LowEnergyPolarizedCompton : public  G4VDiscreteProcess
{  
public:  
  
  G4LowEnergyPolarizedCompton(const G4String& processName = "polarLowEnCompt");
  
  ~G4LowEnergyPolarizedCompton();

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

  G4ThreeVector SetNewPolarization(G4double epsilon, G4double sinSqrTheta, 
				   G4double phi, G4double cosTheta);
  
  G4double SetPhi(G4double, G4double);
  
  void SystemOfRefChange(G4ThreeVector& direction0, G4ThreeVector& direction1, 
			 G4ThreeVector& polarization0, G4ThreeVector& polarization1);

  // hide assignment operator as private 
  G4LowEnergyPolarizedCompton& operator=(const G4LowEnergyPolarizedCompton &right);
  G4LowEnergyPolarizedCompton(const G4LowEnergyPolarizedCompton& );


};

#endif
 








