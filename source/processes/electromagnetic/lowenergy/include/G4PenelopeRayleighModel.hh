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
// $Id: G4PenelopeRayleighModel.hh,v 1.3 2010-04-15 10:02:25 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 13 Oct 2008   L. Pandola   1st implementation. Migration from EM process 
//                            to EM model
// 18 Dec 2009   L. Pandola   Added a dummy ComputeCrossSectioPerAtom() method issueing a
//                            warning if users try to access atomic cross sections via 
//                            G4EmCalculator
//
// -------------------------------------------------------------------
//
// Class description:
// Low Energy Electromagnetic Physics, Rayleigh Scattering
// with Penelope v2001 Model
// -------------------------------------------------------------------

#ifndef G4PENELOPERAYLEIGHMODEL_HH
#define G4PENELOPERAYLEIGHMODEL_HH 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"

class G4ParticleDefinition;
class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4Material;

class G4PenelopeRayleighModel : public G4VEmModel 
{

public:  
  G4PenelopeRayleighModel(const G4ParticleDefinition* p=0,
			  const G4String& processName ="PenRayleigh");
  
  virtual ~G4PenelopeRayleighModel();
  
  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);
  
  virtual G4double CrossSectionPerVolume(const G4Material*,
                                         const G4ParticleDefinition*,
                                         G4double kineticEnergy,
                                         G4double cutEnergy = 0.0,
                                         G4double maxEnergy = DBL_MAX);
  
  //*This is a dummy method. Never inkoved by the tracking, it just issues 
  //*a warning if one tries to get Cross Sections per Atom via the 
  //*G4EmCalculator.
  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                              G4double,
                                              G4double,
                                              G4double,
                                              G4double,
                                              G4double);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);
    
  void SetVerbosityLevel(G4int lev){verboseLevel = lev;};
  G4int GetVerbosityLevel(){return verboseLevel;};
  
protected:
  G4ParticleChangeForGamma* fParticleChange;
  
private:
  G4PenelopeRayleighModel & operator=(const G4PenelopeRayleighModel &right);
  G4PenelopeRayleighModel(const G4PenelopeRayleighModel&);
  
  //Method to initialize sampling 
  void InitialiseSampling();
  std::map <const G4Material*,G4DataVector*> SamplingTable;
  G4DataVector* samplingFunction_x;
  G4DataVector* samplingFunction_xNoLog;

  //Parameters for building the sampling tables
  G4int nPoints; 
  G4double Xhigh;
  G4double Xlow;

  void PrepareConstants();
  
  //Parameters that must be in common between methods
  const G4Material* theMaterial;
  
  //Cross section calculation
  G4double MolecularFormFactor(G4double x);
  G4double DifferentialCrossSection(G4double cosTheta);

  // Energy dependent factor in the differential cross section
  G4double factorE;
  
  //Intrinsic energy limits of the model: cannot be extended by the parent process
  G4double fIntrinsicLowEnergyLimit;
  G4double fIntrinsicHighEnergyLimit;
  
  G4int verboseLevel;
  G4bool isInitialised;
};


#endif

