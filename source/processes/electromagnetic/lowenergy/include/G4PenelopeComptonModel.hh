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
// $Id: G4PenelopeComptonModel.hh,v 1.3 2009-10-21 14:56:47 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 02 Oct 2008   L. Pandola   1st implementation. Migration from EM process 
//                            to EM model
// 21 Oct 2009   L. Pandola   Remove un-necessary methods and variables to handle
//                            AtomicDeexcitationFlag - now demanded to G4VEmModel
//			      Add ActivateAuger() method and deexcitationManager
//
// -------------------------------------------------------------------
//
// Class description:
// Low Energy Electromagnetic Physics, Compton Scattering
// with Penelope Model
// -------------------------------------------------------------------

#ifndef G4PENELOPECOMPTONMODEL_HH
#define G4PENELOPECOMPTONMODEL_HH 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4AtomicDeexcitation.hh"

class G4ParticleDefinition;
class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4Material;

class G4PenelopeComptonModel : public G4VEmModel 
{

public:
  
  G4PenelopeComptonModel(const G4ParticleDefinition* p=0,
			 const G4String& processName ="PenCompton");
  
  virtual ~G4PenelopeComptonModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double ComputeCrossSectionPerAtom(
					      const G4ParticleDefinition*,
					      G4double kinEnergy,
					      G4double Z,
					      G4double A=0,
					      G4double cut=0,
					      G4double emax=DBL_MAX);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

  void SetVerbosityLevel(G4int lev){verboseLevel = lev;};
  G4int GetVerbosityLevel(){return verboseLevel;};

  void ActivateAuger(G4bool);

protected:
  G4ParticleChangeForGamma* fParticleChange;

private:
  //Reads data of atomic shells from database file
  void ReadData();

  //Differential cross section which is numerically integrated
  G4double DifferentialCrossSection (G4double cdt);

  G4PenelopeComptonModel & operator=(const G4PenelopeComptonModel &right);
  G4PenelopeComptonModel(const G4PenelopeComptonModel&);

  //Parameter for the numerical integration of analytical cross section
  G4double energyForIntegration; 
  G4int ZForIntegration;

  //Parameters of atomic shells
  std::map<G4int,G4DataVector*> *ionizationEnergy;
  std::map<G4int,G4DataVector*> *hartreeFunction;
  std::map<G4int,G4DataVector*> *occupationNumber;

  //Intrinsic energy limits of the model: cannot be extended by the parent process
  G4double fIntrinsicLowEnergyLimit;
  G4double fIntrinsicHighEnergyLimit;

  G4int verboseLevel;

  G4bool isInitialised;

  G4AtomicDeexcitation deexcitationManager;

};

#endif

