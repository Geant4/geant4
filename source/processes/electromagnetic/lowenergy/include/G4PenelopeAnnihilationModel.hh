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
// $Id: G4PenelopeAnnihilationModel.hh 74452 2013-10-07 15:08:00Z gcosmo $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 29 Oct 2008   L. Pandola   1st implementation. Migration from EM process 
//                            to EM model
// 02 Oct 2013   L. Pandola   Migration to MT paradigm
//
// -------------------------------------------------------------------
//
// Class description:
// Low Energy Electromagnetic Physics, Positron Annihilation 
// with Penelope Model
// -------------------------------------------------------------------

#ifndef G4PENELOPEANNIHILATIONMODEL_HH
#define G4PENELOPEANNIHILATIONMODEL_HH 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"

class G4ParticleDefinition;
class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4Material;

class G4PenelopeAnnihilationModel : public G4VEmModel 
{
public:
  
  G4PenelopeAnnihilationModel(const G4ParticleDefinition* p=0,
			 const G4String& processName ="PenAnnih");
  
  virtual ~G4PenelopeAnnihilationModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);
  virtual void InitialiseLocal(const G4ParticleDefinition*,
                               G4VEmModel* );

  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
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

protected:
  G4ParticleChangeForGamma* fParticleChange;
  const G4ParticleDefinition* fParticle;

private:
  G4double ComputeCrossSectionPerElectron(G4double energy);

  G4PenelopeAnnihilationModel & operator=(const G4PenelopeAnnihilationModel &right);
  G4PenelopeAnnihilationModel(const G4PenelopeAnnihilationModel&);

  void SetParticle(const G4ParticleDefinition*);

  G4int verboseLevel;
  G4bool isInitialised;

  G4double fIntrinsicLowEnergyLimit;
  G4double fIntrinsicHighEnergyLimit;

  //Stored here so it is not recalculated every time
  static G4double fPielr2;
};

#endif

