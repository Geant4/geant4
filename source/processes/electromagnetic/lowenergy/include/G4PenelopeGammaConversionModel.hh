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
// $Id: G4PenelopeGammaConversionModel.hh,v 1.1 2008-10-28 08:50:21 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 06 Oct 2008   L. Pandola   1st implementation. Migration from EM process 
//                            to EM model
//
// -------------------------------------------------------------------
//
// Class description:
// Low Energy Electromagnetic Physics, Gamma Conversion 
// with Penelope Model
// -------------------------------------------------------------------

#ifndef G4PENELOPEGAMMACONVERSIONMODEL_HH
#define G4PENELOPEGAMMACONVERSIONMODEL_HH 1

#include "globals.hh"
#include "G4VEmModel.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"

class G4ParticleDefinition;
class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4Material;
class G4VCrossSectionHandler;

class G4PenelopeGammaConversionModel : public G4VEmModel 
{

public:
  
  G4PenelopeGammaConversionModel(const G4ParticleDefinition* p=0,
			 const G4String& processName ="PenConversion");
  
  virtual ~G4PenelopeGammaConversionModel();

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

protected:
  G4ParticleChangeForGamma* fParticleChange;

private:
  G4PenelopeGammaConversionModel & operator=(const G4PenelopeGammaConversionModel &right);
  G4PenelopeGammaConversionModel(const G4PenelopeGammaConversionModel&);

  G4double GetScreeningRadius(G4double Z);
  std::vector<G4double>  ScreenFunction(G4double screenVariable);
  G4double CoulombCorrection(G4double ZAlpha);
  G4double LowEnergyCorrection(G4double ZAlpha,G4double eki);

  std::map<G4int,G4double>* fTheScreeningRadii;


  //Intrinsic energy limits of the model: cannot be extended by the parent process
  G4double fIntrinsicLowEnergyLimit;
  G4double fIntrinsicHighEnergyLimit;

  //Use a quicker sampling algorithm if E < smallEnergy
  G4double fSmallEnergy; 

  G4VCrossSectionHandler* crossSectionHandler;

  G4int verboseLevel;
  G4bool isInitialised;
};

#endif

