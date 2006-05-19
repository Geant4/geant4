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
// $Id: G4SynchrotronRadiationInMat.hh,v 1.1 2006-05-19 10:05:28 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      
//      History: 
//      21-5-98  1 version , V. Grichine
//      28-05-01, V.Ivanchenko minor changes to provide ANSI -wall compilation
//      19-05-06, V.Ivanchenko rename from G4SynchrotronRadiation
//
//
// ------------------------------------------------------------

#ifndef G4SynchrotronRadiationInMat_h
#define G4SynchrotronRadiationInMat_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VDiscreteProcess.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4ThreeVector.hh"
#include "G4PropagatorInField.hh"

#include "G4Track.hh"
#include "G4Step.hh"


#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"


#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"


class G4SynchrotronRadiationInMat : public G4VDiscreteProcess
{
public:

  G4SynchrotronRadiationInMat(const G4String& processName =
                              "SynchrotronRadiation",
                              G4ProcessType type = fElectromagnetic);

  virtual ~G4SynchrotronRadiationInMat();

private:

  G4SynchrotronRadiationInMat & operator=(const G4SynchrotronRadiationInMat &right);

  G4SynchrotronRadiationInMat(const G4SynchrotronRadiationInMat&);

public:  /////////////////    Post Step functions  //////////////////////////

  G4double GetMeanFreePath( const G4Track& track,
                            G4double previousStepSize,
                            G4ForceCondition* condition );

  G4VParticleChange *PostStepDoIt( const G4Track& track,
                                      const G4Step& Step    );

  G4double GetPhotonEnergy( const G4Track& trackData,
                               const G4Step&  stepData      );

  G4double GetRandomEnergySR( G4double, G4double );

  G4double GetProbSpectrumSRforInt( G4double );
  G4double GetIntProbSR( G4double );

  G4double GetProbSpectrumSRforEnergy( G4double );
  G4double GetEnergyProbSR( G4double );

  G4double GetIntegrandForAngleK( G4double );
  G4double GetAngleK( G4double );
  G4double GetAngleNumberAtGammaKsi( G4double );



  G4bool IsApplicable(const G4ParticleDefinition&);

  static G4double GetLambdaConst(){ return fLambdaConst; };
  static G4double GetEnergyConst(){ return fEnergyConst; };

  void SetRootNumber(G4int rn){ fRootNumber = rn; };
  void SetVerboseLevel(G4int v){ fVerboseLevel = v; };
  void SetKsi(G4double ksi){ fKsi = ksi; };
  void SetEta(G4double eta){ fEta = eta; };
  void SetPsiGamma(G4double psg){ fPsiGamma = psg; };
  void SetOrderAngleK(G4double ord){ fOrderAngleK = ord; }; // should be 1/3 or 2/3

  protected:

  private:

  static const G4double fLambdaConst;

  static const G4double fEnergyConst;

  static const G4double fIntegralProbabilityOfSR[200];


  const G4double
  LowestKineticEnergy;   // low  energy limit of the cross-section formula

  const G4double
  HighestKineticEnergy;  // high energy limit of the cross-section formula

  G4int TotBin;          // number of bins in the tables

  G4double CutInRange;

  const G4ParticleDefinition* theGamma;
  const G4ParticleDefinition* theElectron;
  const G4ParticleDefinition* thePositron;

  const G4double* GammaCutInKineticEnergy;
  const G4double* ElectronCutInKineticEnergy;
  const G4double* PositronCutInKineticEnergy;
  const G4double* ParticleCutInKineticEnergy;


  G4double GammaCutInKineticEnergyNow;
  G4double ElectronCutInKineticEnergyNow;
  G4double PositronCutInKineticEnergyNow;
  G4double ParticleCutInKineticEnergyNow;

  G4double fAlpha;
  G4int    fRootNumber;
  G4double fKsi;             // omega/omega_c
  G4double fPsiGamma;        // Psi-angle*gamma
  G4double fEta;             //
  G4double fOrderAngleK;     // 1/3 or 2/3


  G4int    fVerboseLevel;
  G4PropagatorInField* fFieldPropagator;

};

//////////////////////////  INLINE METHODS  /////////////////////////////

inline G4bool
G4SynchrotronRadiationInMat::IsApplicable( const G4ParticleDefinition& particle )
{

  return ( ( &particle == (const G4ParticleDefinition *)theElectron ) ||
            ( &particle == (const G4ParticleDefinition *)thePositron )    );

  // return ( particle.GetPDGCharge() != 0.0 );
}

#endif  // end of G4SynchrotronRadiationInMat.hh

