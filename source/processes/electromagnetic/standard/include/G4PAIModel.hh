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
// $Id: G4PAIModel.hh,v 1.22 2009-02-19 19:17:50 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4PAIModel
//
// Author:        V. Grichine based on Vladimir Ivanchenko  code
//
// Creation date: 05.10.2003
//
// Modifications:
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 26-09-07 Fixed tmax computation (V.Ivantchenko)
//
//
// Class Description:
//
// Implementation of PAI model of energy loss and
// delta-electron production by heavy charged particles

// -------------------------------------------------------------------
//

#ifndef G4PAIModel_h
#define G4PAIModel_h 1

#include <vector>
#include "G4VEmModel.hh"
#include "globals.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4PAIySection.hh"

class G4PhysicsLogVector;
class G4PhysicsTable;
class G4Region;
class G4MaterialCutsCouple;
class G4ParticleChangeForLoss;

class G4PAIModel : public G4VEmModel, public G4VEmFluctuationModel
{

public:

  G4PAIModel(const G4ParticleDefinition* p = 0, const G4String& nam = "PAI");

  virtual ~G4PAIModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual void InitialiseMe(const G4ParticleDefinition*);

  virtual G4double ComputeDEDXPerVolume(const G4Material*,
			       const G4ParticleDefinition*,
			       G4double kineticEnergy,
			       G4double cutEnergy);

  virtual G4double CrossSectionPerVolume(const G4Material*,
				const G4ParticleDefinition*,
				G4double kineticEnergy,
				G4double cutEnergy,
				G4double maxEnergy);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

  virtual G4double SampleFluctuations(const G4Material*,
				      const G4DynamicParticle*,
				      G4double&,
				      G4double&,
				      G4double&);

  virtual G4double Dispersion(    const G4Material*,
				  const G4DynamicParticle*,
				  G4double&,
				  G4double&);

  void     DefineForRegion(const G4Region* r) ;
  void     ComputeSandiaPhotoAbsCof();
  void     BuildPAIonisationTable();
  void     BuildLambdaVector();

  G4double GetdNdxCut( G4int iPlace, G4double transferCut);
  G4double GetdEdxCut( G4int iPlace, G4double transferCut);
  G4double GetPostStepTransfer( G4double scaledTkin );
  G4double GetEnergyTransfer( G4int iPlace,
			      G4double position,
			      G4int iTransfer );

  void SetVerboseLevel(G4int verbose){fVerbose=verbose;};

protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
                                    G4double kinEnergy);

private:

  void SetParticle(const G4ParticleDefinition* p);

  // hide assignment operator 
  G4PAIModel & operator=(const  G4PAIModel &right);
  G4PAIModel(const  G4PAIModel&);

  // The vector over proton kinetic energies: the range of gammas
  G4int                fVerbose; 
  G4double             fLowestGamma;
  G4double             fHighestGamma;
  G4double             fLowestKineticEnergy;
  G4double             fHighestKineticEnergy;
  G4int                fTotBin;
  G4int                fMeanNumber;
  G4PhysicsLogVector*  fParticleEnergyVector ;
  G4PAIySection        fPAIySection;

  // vectors

  G4PhysicsTable*                    fPAItransferTable;
  std::vector<G4PhysicsTable*>       fPAIxscBank;

  G4PhysicsTable*                    fPAIdEdxTable;
  std::vector<G4PhysicsTable*>       fPAIdEdxBank;

  std::vector<const G4MaterialCutsCouple*> fMaterialCutsCoupleVector;
  std::vector<const G4Region*>       fPAIRegionVector;

  const G4MaterialCutsCouple*        fCutCouple;
  const G4Material*                  fMaterial;
  G4double                           fDeltaCutInKinEnergy; 

  size_t                             fMatIndex ;  
  G4double**                         fSandiaPhotoAbsCof ;
  G4int                              fSandiaIntervalNumber ;

  G4PhysicsLogVector*                fdEdxVector ;
  std::vector<G4PhysicsLogVector*>   fdEdxTable ;

  G4PhysicsLogVector*                fLambdaVector ;
  std::vector<G4PhysicsLogVector*>   fLambdaTable ;

  G4PhysicsLogVector*                fdNdxCutVector ;
  std::vector<G4PhysicsLogVector*>   fdNdxCutTable ;

  const G4ParticleDefinition* fParticle;
  const G4ParticleDefinition* fElectron;
  const G4ParticleDefinition* fPositron;
  G4ParticleChangeForLoss*    fParticleChange;

  G4double fMass;
  G4double fSpin;
  G4double fChargeSquare;
  G4double fRatio;
  G4double fHighKinEnergy;
  G4double fLowKinEnergy;
  G4double fTwoln10;
  G4double fBg2lim; 
  G4double fTaulim;
  G4double fQc;

  G4bool   isInitialised;
};

#endif







