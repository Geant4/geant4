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
// $Id: G4PenelopeIonisationModel.hh,v 1.5 2010-04-15 10:02:10 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 26 Nov 2008   L. Pandola   1st implementation. Migration from EM process 
//                            to EM model. Physics is unchanged.
// 21 Oct 2009   L. Pandola   Remove un-necessary methods and variables to handle 
//                            AtomicDeexcitationFlag - now demanded to G4VEmModel
//			      Add ActivateAuger() method
// 29 Mar 2010   L. Pandola   Added a dummy ComputeCrossSectioPerAtom() method issueing a
//                            warning if users try to access atomic cross sections via 
//                            G4EmCalculator
// 15 Apr 2010   L. Pandola   Implemented model's own version of MinEnergyCut()
//
// -------------------------------------------------------------------
//
// Class description:
// Low Energy Electromagnetic Physics, e+ and e- ionisation
// with Penelope Model
// -------------------------------------------------------------------

#ifndef G4PENELOPEIONISATIONMODEL_HH
#define G4PENELOPEIONISATIONMODEL_HH 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4VCrossSectionHandler.hh"
#include "G4PhysicsLogVector.hh"
#include "G4AtomicDeexcitation.hh"

class G4ParticleDefinition;
class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4Material;
class G4VEMDataSet;

class G4PenelopeIonisationModel : public G4VEmModel 
{

public:
  
  G4PenelopeIonisationModel(const G4ParticleDefinition* p=0,
			 const G4String& processName ="PenelopeIoni");
  
  virtual ~G4PenelopeIonisationModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

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
                                         const G4ParticleDefinition* theParticle,
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

  void ActivateAuger(G4bool);

protected:
  G4ParticleChangeForLoss* fParticleChange;

private:
 
  G4PenelopeIonisationModel & operator=(const G4PenelopeIonisationModel &right);
  G4PenelopeIonisationModel(const G4PenelopeIonisationModel&);

 
  //Intrinsic energy limits of the model: cannot be extended by the parent process
  G4double fIntrinsicLowEnergyLimit;
  G4double fIntrinsicHighEnergyLimit;

  G4int verboseLevel;

  G4bool isInitialised;
 
  G4double CalculateDeltaFermi(G4double kinEnergy ,G4int Z,
			       G4double electronVolumeDensity);
  	
  //Methods and variables to calculate final state
  void CalculateDiscreteForElectrons(G4double kinEnergy,G4double cutoffEnergy,
				     G4int Z,G4double electronVolumeDensity);
  void CalculateDiscreteForPositrons(G4double kinEnergy,G4double cutoffEnergy,
			     G4int Z,G4double electronVolumeDensity);

  G4AtomicDeexcitation deexcitationManager;
  G4double kineticEnergy1;
  G4double cosThetaPrimary;
  G4double energySecondary;
  G4double cosThetaSecondary;
  G4int iOsc;				   

  //These methods are used to calculate the hard-cross section (namely they 
  //return the hard/total cross section)
  G4double CalculateCrossSectionsRatio(G4double kinEnergy,
				       G4double cutoffEnergy,
				       G4int Z, 
				       G4double electronVolumeDensity,
				       const G4ParticleDefinition*);
  //In fact the total cross section (hard+soft) is read from file
  //The following methods give the cross section contribution (hard and soft) from each 
  //individual oscillator
  std::pair<G4double,G4double> CrossSectionsRatioForElectrons(G4double kineticEnergy,
							      G4double resEnergy,
							      G4double densityCorrection,
							      G4double cutoffEnergy);

  std::pair<G4double,G4double> CrossSectionsRatioForPositrons(G4double kineticEnergy,
							      G4double resEnergy,
							      G4double densityCorrection,
							      G4double cutoffEnergy);
  
  G4VCrossSectionHandler* crossSectionHandler;
  
  //These methods are used to calculate the stopping power up to the cutoff
  //for each individual oscillator
  G4double ComputeStoppingPowerForElectrons(G4double kinEnergy,
					    G4double cutEnergy,
					    G4double deltaFermi,
					    G4double resEnergy);

  G4double ComputeStoppingPowerForPositrons(G4double kinEnergy,
					    G4double cutEnergy,
					    G4double deltaFermi,
					    G4double resEnergy);
  
  
  //Parameters of atomic shells
  void ReadData();
  std::map<G4int,G4DataVector*> *ionizationEnergy;
  std::map<G4int,G4DataVector*> *resonanceEnergy;
  std::map<G4int,G4DataVector*> *occupationNumber;
  std::map<G4int,G4DataVector*> *shellFlag;
  
  //Mean free path table. This will become obsolete! For now I need something to store 
  //cross sections and to sample a random atom
  std::vector<G4VEMDataSet*>* theXSTable;
  std::vector<G4VEMDataSet*>* BuildCrossSectionTable(const G4ParticleDefinition*);
  G4int SampleRandomAtom(const G4MaterialCutsCouple*,G4double energy) const;

};

#endif

