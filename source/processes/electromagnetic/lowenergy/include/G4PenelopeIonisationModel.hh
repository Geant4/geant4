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
// $Id: G4PenelopeIonisationModel.hh 91228 2015-06-26 10:37:01Z gcosmo $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 30 Mar 2010   L. Pandola   1st implementation. 
// 25 May 2011   L. Pandola   Renamed (make v2008 as default Penelope)
// 09 Mar 2012   L. Pandola   Moved the management and calculation of 
//                            cross sections to a separate class
// 07 Oct 2013   L. Pandola   Migration to MT
// 23 Jun 2015   L. Pandola   Added private member to store the PIXE flag
//
// -------------------------------------------------------------------
//
// Class description:
// Low Energy Electromagnetic Physics, e+ and e- ionisation
// with Penelope Model, version 2008
// -------------------------------------------------------------------

#ifndef G4PENELOPEIONISATIONMODEL_HH
#define G4PENELOPEIONISATIONMODEL_HH 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4VAtomDeexcitation.hh"

class G4PhysicsFreeVector;
class G4PhysicsLogVector;
class G4ParticleDefinition;
class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4Material;
class G4PenelopeOscillatorManager;
class G4PenelopeOscillator;
class G4PenelopeCrossSection;
class G4PenelopeIonisationXSHandler;

class G4PenelopeIonisationModel : public G4VEmModel 
{

public:
  
  G4PenelopeIonisationModel(const G4ParticleDefinition* p=0,
			    const G4String& processName ="PenIoni");
  
  virtual ~G4PenelopeIonisationModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);
  virtual void InitialiseLocal(const G4ParticleDefinition*,
			       G4VEmModel*);

  //*This is a dummy method. Never inkoved by the tracking, it just issues 
  //*a warning if one tries to get Cross Sections per Atom via the 
  //*G4EmCalculator.
  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                              G4double,
                                              G4double,
                                              G4double,
                                              G4double,
                                              G4double);

  virtual G4double CrossSectionPerVolume(const G4Material* material,
                                         const G4ParticleDefinition* 
					 theParticle,
                                         G4double kineticEnergy,
                                         G4double cutEnergy,
                                         G4double maxEnergy = DBL_MAX);
					 
  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);
				   
  virtual G4double ComputeDEDXPerVolume(const G4Material*,
                               const G4ParticleDefinition*,
                               G4double kineticEnergy,
                               G4double cutEnergy);

  // Min cut in kinetic energy allowed by the model
  virtual G4double MinEnergyCut(const G4ParticleDefinition*,
                                const G4MaterialCutsCouple*);		

  void SetVerbosityLevel(G4int lev){verboseLevel = lev;};
  G4int GetVerbosityLevel(){return verboseLevel;};

protected:
  G4ParticleChangeForLoss* fParticleChange;
  const G4ParticleDefinition* fParticle;

private:
 
  G4PenelopeIonisationModel & operator=(const G4PenelopeIonisationModel &right);
  G4PenelopeIonisationModel(const G4PenelopeIonisationModel&);

  void SetParticle(const G4ParticleDefinition*);

  void SampleFinalStateElectron(const G4Material*,
				G4double cutEnergy,
				G4double kineticEnergy);
  void SampleFinalStatePositron(const G4Material*,
				G4double cutEnergy,
				G4double kineticEnergy);

  //Intrinsic energy limits of the model: cannot be extended by the parent process
  G4double fIntrinsicLowEnergyLimit;
  G4double fIntrinsicHighEnergyLimit;

  G4int verboseLevel;

  G4bool isInitialised;
  G4VAtomDeexcitation* fAtomDeexcitation;
  G4bool fPIXEflag;


  G4double kineticEnergy1;
  G4double cosThetaPrimary;
  G4double energySecondary;
  G4double cosThetaSecondary;
  G4int targetOscillator;				   

  G4PenelopeOscillatorManager* oscManager;
  G4PenelopeIonisationXSHandler* theCrossSectionHandler;

  size_t nBins;

  //Used only for G4EmCalculator and Unit Tests
  G4bool fLocalTable;

};

#endif

