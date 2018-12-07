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
// -------------------------------------------------------------------
//
// Author: A. Forti
//
// History:
// -----------
// 02 Mar 1999  A. Forti        1st implementation
// 27 Sep 2001  V. Ivanchenko   Major revision according to a design iteration
// 10 Oct 2001  M.G. Pia        Revision to improve code quality 
//                              and consistency with design
// 29 Nov 2001  V.Ivanchenko    New parametrisation of EEDL data
// 21 Feb 2003  V.Ivanchenko    Energy bins for spectrum are defined here
// 24 Mar 2003  P.Rodrigues     Changes to accommodate new angular generators
// 06 Nov 2003  MGP             Overloaded SetAngularGenerator
//
// -------------------------------------------------------------------

// Class description:
// Low Energy electromagnetic process, electron Bremsstrahlung
// based on the data of the EEDL database. Details are described
// in the Physics Reference Manual.
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// --------------------------------------------------------------


#ifndef G4RDLOWENERGYBREMSSTRAHLUNG_HH
#define G4RDLOWENERGYBREMSSTRAHLUNG_HH 1

#include "G4eLowEnergyLoss.hh"
#include "G4DataVector.hh"

class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4VParticleChange;
class G4RDVEMDataSet;
class G4RDVEnergySpectrum;
class G4RDVCrossSectionHandler;
class G4RDVBremAngularDistribution;

class G4LowEnergyBremsstrahlung : public G4eLowEnergyLoss
{ 
public:
 
  G4LowEnergyBremsstrahlung(const G4String& processName = "LowEnBrem");
  
  
  //  G4LowEnergyBremsstrahlung(const G4String& processName = "LowEnBrem",
  //			    G4RDVBremAngularDistribution* distribution = 0);
  
  ~G4LowEnergyBremsstrahlung();
  
  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void BuildPhysicsTable(const G4ParticleDefinition& particleType);
  
  G4VParticleChange* PostStepDoIt(const G4Track& track,         
				  const G4Step& step);
  
  void SetCutForLowEnSecPhotons(G4double cut);
  
  void SetAngularGenerator(G4RDVBremAngularDistribution* distribution);

  void SetAngularGenerator(const G4String& name);
  
  void PrintInfoDefinition();
        

protected:

  G4double GetMeanFreePath(const G4Track& track,
			   G4double previousStepSize,
			   G4ForceCondition* condition );
 
private:

  // Hide copy constructor and assignment operator as private 
  G4LowEnergyBremsstrahlung(const G4LowEnergyBremsstrahlung& );
  G4LowEnergyBremsstrahlung& operator = (const G4LowEnergyBremsstrahlung& right);
  
  void BuildLossTable(const G4ParticleDefinition& ParticleType);
  
  G4RDVCrossSectionHandler* crossSectionHandler;
  G4RDVEMDataSet* theMeanFreePath;
  G4RDVEnergySpectrum* energySpectrum;
  G4DataVector  energyBins;
  G4RDVBremAngularDistribution* angularDistribution;
  G4RDVBremAngularDistribution* TsaiAngularDistribution;
  G4String generatorName;

  // Lower limit for generation of gamma in this model
  G4DataVector cutForSecondaryPhotons;
  G4double cutForPhotons;
};



#endif











