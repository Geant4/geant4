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
// $Id: G4EmParameters.hh 66885 2013-01-16 17:37:13Z gunter $
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4EmParameters
//
// Author:        Vladimir Ivanchenko for migration to MT
//                  
//
// Creation date: 17.05.2013
//
// Modifications:
//
//
// Class Description:
//
// A utility static class, responsable for keeping parameters
// for all EM physics processes and models.
//
// It is initialized by the master thread but can be updated 
// at any moment. Parameters may be used in run time or at 
// initialisation
//
// -------------------------------------------------------------------
//

#ifndef G4EmParameters_h
#define G4EmParameters_h 1

#include "globals.hh"
#include "G4MscStepLimitType.hh"

class G4EmParametersMessenger;

class G4EmParameters
{
public:

  static G4EmParameters* Instance();

  ~G4EmParameters();

  void SetLossFluctuations(G4bool val);
  G4bool LossFluctuation() const;

  void SetBuildCSDARange(G4bool val);
  G4bool BuildCSDARange() const;

  void SetLPM(G4bool val);
  G4bool LPM() const;

  void SetSpline(G4bool val);
  G4bool Spline() const;

  void SetUseCutAsFinalRange(G4bool val);
  G4bool UseCutAsFinalRange() const;

  void SetApplyCuts(G4bool val);
  G4bool ApplyCuts() const;

  void SetFluo(G4bool val);
  G4bool Fluo() const;

  void SetAuger(G4bool val);
  G4bool Auger() const;

  void SetPixe(G4bool val);
  G4bool Pixe() const;

  void SetDeexcitationIgnoreCut(G4bool val);
  G4bool DeexcitationIgnoreCut() const;

  void SetLateralDisplacement(G4bool val);
  G4bool LateralDisplacement() const;

  void SetMuHadLateralDisplacement(G4bool val);
  G4bool MuHadLateralDisplacement() const;

  void SetLatDisplacementBeyondSafety(G4bool val);
  G4bool LatDisplacementBeyondSafety() const;

  void ActivateAngularGeneratorForIonisation(G4bool val);
  G4bool UseAngularGeneratorForIonisation() const;

  void SetMinSubRange(G4double val);
  G4double MinSubRange() const;

  void SetMinEnergy(G4double val);
  G4double MinKinEnergy() const;

  void SetMaxEnergy(G4double val);
  G4double MaxKinEnergy() const;

  void SetMaxEnergyForCSDARange(G4double val);
  G4double MaxEnergyForCSDARange() const;

  void SetLinearLossLimit(G4double val);
  G4double LinearLossLimit() const;

  void SetBremsstrahlungTh(G4double val);
  G4double BremsstrahlungTh() const;

  void SetLambdaFactor(G4double val);
  G4double LambdaFactor() const;

  void SetFactorForAngleLimit(G4double val);
  G4double FactorForAngleLimit() const;

  void SetMscThetaLimit(G4double val);
  G4double MscThetaLimit() const;

  void SetMscRangeFactor(G4double val);
  G4double MscRangeFactor() const;

  void SetMscGeomFactor(G4double val);
  G4double MscGeomFactor() const;

  void SetMscSkin(G4double val);
  G4double MscSkin() const;

  void SetNumberOfBins(G4int val);
  G4int NumberOfBins() const;

  void SetNumberOfBinsPerDecade(G4int val);
  G4int NumberOfBinsPerDecade() const;

  void SetVerbose(G4int val);
  G4int Verbose() const;

  void SetWorkerVerbose(G4int val);
  G4int WorkerVerbose() const;

  void SetMscStepLimitType(G4MscStepLimitType val);
  G4MscStepLimitType MscStepLimitType() const;

private:

  G4EmParameters();

  void PrintWarning(G4ExceptionDescription& ed);

  static G4EmParameters* theInstance;

  G4EmParametersMessenger* theMessenger;

  G4bool lossFluctuation;
  G4bool buildCSDARange;
  G4bool flagLPM;
  G4bool spline;
  G4bool finalRange;
  G4bool applyCuts;
  G4bool fluo;
  G4bool auger;
  G4bool pixe;
  G4bool deexIgnoreCut;
  G4bool lateralDisplacement;
  G4bool muhadLateralDisplacement;
  G4bool latDisplacementBeyondSafety;
  G4bool useAngGeneratorForIonisation;

  G4double minSubRange;
  G4double minKinEnergy;
  G4double maxKinEnergy;
  G4double maxKinEnergyCSDA;
  G4double linLossLimit;
  G4double bremsTh;
  G4double lambdaFactor;
  G4double factorForAngleLimit;
  G4double thetaLimit;
  G4double rangeFactor;
  G4double geomFactor;
  G4double skin;

  G4int nbins;
  G4int nbinsPerDecade;
  G4int verbose;
  G4int workerVerbose;

  G4MscStepLimitType mscStepLimit;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

