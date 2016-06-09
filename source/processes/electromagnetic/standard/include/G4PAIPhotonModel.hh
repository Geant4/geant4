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
// $Id: G4PAIPhotonModel.hh,v 1.12 2009-02-19 19:17:50 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4PAIPhotonModel
//
// Author:        V. Grichine based on Vladimir Ivanchenko  code
//
// Creation date: 05.10.2003
//
// Modifications:
// 11.04.05 Major optimisation of internal interfaces (V.Ivantchenko)
// 26-09-07 Fixed tmax computation (V.Ivantchenko)
//
//
// Class Description:
//
// Implementation of PAI model of energy loss and
// delta-electron production by heavy charged particles

// -------------------------------------------------------------------
//

#ifndef G4PAIPhotonModel_h
#define G4PAIPhotonModel_h 1

#include <vector>
#include "G4VEmModel.hh"
#include "globals.hh"
#include "G4VEmFluctuationModel.hh"

class G4PhysicsLogVector;
class G4PhysicsTable;
class G4Region;
class G4MaterialCutsCouple;
class G4ParticleChangeForLoss;

class G4PAIPhotonModel : public G4VEmModel, public G4VEmFluctuationModel
{

public:

  G4PAIPhotonModel(const G4ParticleDefinition* p = 0, const G4String& nam = "PAIPhoton");

  virtual ~G4PAIPhotonModel();

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
  void     BuildLambdaVector(const G4MaterialCutsCouple* matCutsCouple);

  G4double GetdNdxCut( G4int iPlace, G4double transferCut);
  G4double GetdNdxPhotonCut( G4int iPlace, G4double transferCut);
  G4double GetdNdxPlasmonCut( G4int iPlace, G4double transferCut);

  G4double GetdEdxCut( G4int iPlace, G4double transferCut);

  G4double GetPostStepTransfer(G4PhysicsTable*, G4PhysicsLogVector*,
                               G4int iPlace, G4double scaledTkin );
  G4double GetAlongStepTransfer(G4PhysicsTable*, G4PhysicsLogVector*,
                               G4int iPlace, G4double scaledTkin,G4double step, G4double cof );
  G4double GetEnergyTransfer(G4PhysicsTable*, G4int iPlace,
                             G4double position, G4int iTransfer );

protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
			      G4double kinEnergy);

private:

  void SetParticle(const G4ParticleDefinition* p);

  // hide assignment operator 
  G4PAIPhotonModel & operator=(const  G4PAIPhotonModel &right);
  G4PAIPhotonModel(const  G4PAIPhotonModel&);

  // The vector over proton kinetic energies: the range of gammas

  G4double             fLowestKineticEnergy;
  G4double             fHighestKineticEnergy;
  G4int                fTotBin;
  G4int                fMeanNumber;
  G4int                fVerbose; 
  G4PhysicsLogVector*  fProtonEnergyVector ;

  // vectors

  G4PhysicsTable*                    fPAItransferTable;
  std::vector<G4PhysicsTable*>       fPAIxscBank;

  G4PhysicsTable*                    fPAIphotonTable;
  std::vector<G4PhysicsTable*>       fPAIphotonBank;

  G4PhysicsTable*                    fPAIplasmonTable;
  std::vector<G4PhysicsTable*>       fPAIplasmonBank;

  G4PhysicsTable*                    fPAIdEdxTable;
  std::vector<G4PhysicsTable*>       fPAIdEdxBank;

  std::vector<const G4MaterialCutsCouple*> fMaterialCutsCoupleVector;
  std::vector<const G4Region*>       fPAIRegionVector;

  size_t                             fMatIndex ;  
  G4double**                         fSandiaPhotoAbsCof ;
  G4int                              fSandiaIntervalNumber ;

  G4PhysicsLogVector*              fdEdxVector ;
  std::vector<G4PhysicsLogVector*> fdEdxTable ;

  G4PhysicsLogVector*              fLambdaVector ;
  std::vector<G4PhysicsLogVector*> fLambdaTable ;

  G4PhysicsLogVector*              fdNdxCutVector ;
  std::vector<G4PhysicsLogVector*> fdNdxCutTable ;

  G4PhysicsLogVector*              fdNdxCutPhotonVector ;
  std::vector<G4PhysicsLogVector*> fdNdxCutPhotonTable ;

  G4PhysicsLogVector*              fdNdxCutPlasmonVector ;
  std::vector<G4PhysicsLogVector*> fdNdxCutPlasmonTable ;


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







