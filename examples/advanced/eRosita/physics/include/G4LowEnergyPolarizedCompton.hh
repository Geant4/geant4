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
// 17 October 2001 - F.Longo  Major revision according to design iteration
//
// 21 February 2002 - F.Longo Revisions with A.Zoglauer and G.Depaola
//                            - better description of parallelism
//                            - system of ref change method improved
//
//
// Class description:
// Low Energy electromagnetic process, Polarised Compton scattering
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// ------------------------------------------------------------

#ifndef G4RDLOWENERGYPOLARIZEDCOMPTON_H
#define G4RDLOWENERGYPOLARIZEDCOMPTON_H 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"

// Doppler Broadening

#include "G4RDShellData.hh"
#include "G4RDDopplerProfile.hh"


class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4RDVEMDataSet;
class G4RDVCrossSectionHandler;
class G4RDVRangeTest;

class G4LowEnergyPolarizedCompton : public  G4VDiscreteProcess
{  
public:  
  
  G4LowEnergyPolarizedCompton(const G4String& processName = "polarLowEnCompt");
  
  ~G4LowEnergyPolarizedCompton();

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

  G4LowEnergyPolarizedCompton& operator=(const G4LowEnergyPolarizedCompton& 
					 right);
  G4LowEnergyPolarizedCompton(const G4LowEnergyPolarizedCompton& );
  
  G4double lowEnergyLimit;  // low energy limit  applied to the process
  G4double highEnergyLimit; // high energy limit applied to the process
  
  G4RDVEMDataSet* meanFreePathTable;
  G4RDVEMDataSet* scatterFunctionData;

  G4RDVCrossSectionHandler* crossSectionHandler;
  G4RDVRangeTest* rangeTest;

  const G4double intrinsicLowEnergyLimit; // intrinsic validity range
  const G4double intrinsicHighEnergyLimit;

  // specific methods for polarization 
  
  G4ThreeVector GetRandomPolarization(G4ThreeVector& direction0); // Random Polarization
  G4ThreeVector GetPerpendicularPolarization(const G4ThreeVector& direction0, const G4ThreeVector& polarization0) const;
  
  G4ThreeVector SetPerpendicularVector(G4ThreeVector& a); // temporary
  G4ThreeVector SetNewPolarization(G4double epsilon, G4double sinSqrTheta, 
				   G4double phi, G4double cosTheta);
  G4double SetPhi(G4double, G4double);
  
  void SystemOfRefChange(G4ThreeVector& direction0, G4ThreeVector& direction1, 
			 G4ThreeVector& polarization0, G4ThreeVector& polarization1);
  
  // Doppler Broadening
 
  G4RDShellData shellData;
  G4RDDopplerProfile profileData;


};

#endif
 








