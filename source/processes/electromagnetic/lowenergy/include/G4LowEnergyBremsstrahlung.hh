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
// $Id: G4LowEnergyBremsstrahlung.hh,v 1.22 2001-10-10 09:49:29 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: A. Forti
//
// History:
// -----------
// 02 Mar 1999  A. Forti        1st implementation
// 27 Sep 2001  V. Ivanchenko   Major revision according to a design iteration
// 10 Oct 2001  M.G. Pia        Revision to improve code quality and consistency with design
//
// -------------------------------------------------------------------

// Class description:
// Low Energy electromagnetic process, electron Bremsstrahlung
// Further documentation available from http://www.ge.infn.it/geant4/lowE
//
// Class Description: End 

// --------------------------------------------------------------
//

#ifndef G4LOWENERGYBREMSSTRAHLUNG_HH
#define G4LOWENERGYBREMSSTRAHLUNG_HH 1

#include "G4eLowEnergyLoss.hh"
#include "G4DataVector.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4VEMDataSet;
class G4VEnergySpectrum;
class G4VCrossSectionHandler;


class G4LowEnergyBremsstrahlung : public G4eLowEnergyLoss

{ 
public:
 
  G4LowEnergyBremsstrahlung(const G4String& processName = "eLowEnergyBrem");
  
  ~G4LowEnergyBremsstrahlung();
  
  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void PrintInfoDefinition();
  
  void BuildPhysicsTable(const G4ParticleDefinition& particleType);
  
  G4VParticleChange* PostStepDoIt(const G4Track& track,         
				  const G4Step& step);                 
      

protected:

  G4double GetMeanFreePath(const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition );
 
private:

  // Hide copy constructor and assignment operator as private 
  G4LowEnergyBremsstrahlung(const G4LowEnergyBremsstrahlung& );
  G4LowEnergyBremsstrahlung& operator = (const G4LowEnergyBremsstrahlung& right);
  
  void BuildLossTable(const G4ParticleDefinition& ParticleType);
  
  G4VCrossSectionHandler* crossSectionHandler;
  G4VEMDataSet* theMeanFreePath;
  G4VEnergySpectrum* energySpectrum;

  // Lower limit for generation of gamma in this model
  G4DataVector cutForSecondaryPhotons;

};


  
#endif
 










