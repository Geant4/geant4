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
#include "G4ios.hh"
#include "G4MscStepLimitType.hh"
#include "G4NuclearFormfactorType.hh"
#include "G4EmSaturation.hh"
#include "G4Threading.hh"
#include <vector>

class G4EmParametersMessenger;
class G4VEnergyLossProcess;
class G4VEmProcess;
class G4VAtomDeexcitation;
class G4StateManager;

class G4EmParameters
{
public:

  static G4EmParameters* Instance();

  ~G4EmParameters();

  void SetDefaults();

  // printing
  std::ostream& StreamInfo(std::ostream& os) const;
  void Dump() const;
  friend std::ostream& operator<< (std::ostream& os, const G4EmParameters&);

  // boolean flags
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

  void SetBeardenFluoDir(G4bool val);
  G4bool BeardenFluoDir() const;

  void SetAuger(G4bool val);
  G4bool Auger() const;

  void SetAugerCascade(G4bool val);
  G4bool AugerCascade() const;

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

  void SetUseMottCorrection(G4bool val);
  G4bool UseMottCorrection() const;

  void SetIntegral(G4bool val);
  G4bool Integral() const;

  void SetBirksActive(G4bool val);
  G4bool BirksActive() const;

  void SetEmSaturation(G4EmSaturation*);
  G4EmSaturation* GetEmSaturation();

  // double parameters with values
  void SetMinSubRange(G4double val);
  G4double MinSubRange() const;

  void SetMinEnergy(G4double val);
  G4double MinKinEnergy() const;

  void SetMaxEnergy(G4double val);
  G4double MaxKinEnergy() const;

  void SetMaxEnergyForCSDARange(G4double val);
  G4double MaxEnergyForCSDARange() const;

  void SetLowestElectronEnergy(G4double val);
  G4double LowestElectronEnergy() const;

  void SetLowestMuHadEnergy(G4double val);
  G4double LowestMuHadEnergy() const;

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

  void SetMscMuHadRangeFactor(G4double val);
  G4double MscMuHadRangeFactor() const;

  void SetMscGeomFactor(G4double val);
  G4double MscGeomFactor() const;

  void SetMscSkin(G4double val);
  G4double MscSkin() const;

  void SetScreeningFactor(G4double val);
  G4double ScreeningFactor() const;

  void SetStepFunction(G4double v1, G4double v2);

  void SetStepFunctionMuHad(G4double v1, G4double v2);

  // integer parameters 
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

  void SetMscMuHadStepLimitType(G4MscStepLimitType val);
  G4MscStepLimitType MscMuHadStepLimitType() const;

  void SetNuclearFormfactorType(G4NuclearFormfactorType val);
  G4NuclearFormfactorType NuclearFormfactorType() const;

  // string parameters 
  void SetPIXECrossSectionModel(const G4String&);
  const G4String& PIXECrossSectionModel();

  void SetPIXEElectronCrossSectionModel(const G4String&);
  const G4String& PIXEElectronCrossSectionModel();

  // parameters per region or per process 
  void AddPAIModel(const G4String& particle,
                   const G4String& region,
                   const G4String& type);
  const std::vector<G4String>& ParticlesPAI() const;
  const std::vector<G4String>& RegionsPAI() const;
  const std::vector<G4String>& TypesPAI() const;

  void AddMicroElec(const G4String& region);
  const std::vector<G4String>& RegionsMicroElec() const;

  void AddDNA(const G4String& region, const G4String& type);
  const std::vector<G4String>& RegionsDNA() const;
  const std::vector<G4String>& TypesDNA() const;

  void AddMsc(const G4String& region, const G4String& type);
  const std::vector<G4String>& RegionsMsc() const;
  const std::vector<G4String>& TypesMsc() const;

  void AddPhysics(const G4String& region, const G4String& type);
  const std::vector<G4String>& RegionsPhysics() const;
  const std::vector<G4String>& TypesPhysics() const;

  void SetSubCutoff(G4bool val, const G4String& region = "");

  void SetDeexActiveRegion(const G4String& region, G4bool fdeex,
			   G4bool fauger, G4bool fpixe);

  void SetProcessBiasingFactor(const G4String& procname, 
                               G4double val, G4bool wflag);

  void ActivateForcedInteraction(const G4String& procname, 
                                 const G4String& region,
                                 G4double length, 
                                 G4bool wflag);

  void ActivateSecondaryBiasing(const G4String& name,
				const G4String& region, 
				G4double factor,
				G4double energyLimit);

  // initialisation methods
  void DefineRegParamForLoss(G4VEnergyLossProcess*, 
                             G4bool isElectron) const;
  void DefineRegParamForEM(G4VEmProcess*) const;
  void DefineRegParamForDeex(G4VAtomDeexcitation*) const;

private:

  G4EmParameters(G4EmParameters &) = delete;
  G4EmParameters & operator=(const G4EmParameters &right) = delete;  

  G4EmParameters();

  void Initialise();

  G4bool IsLocked() const;

  G4String CheckRegion(const G4String&) const;

  void PrintWarning(G4ExceptionDescription& ed) const;

  static G4EmParameters* theInstance;

  G4EmParametersMessenger* theMessenger;

  G4StateManager* fStateManager;

  G4EmSaturation* emSaturation;

  G4bool lossFluctuation;
  G4bool buildCSDARange;
  G4bool flagLPM;
  G4bool spline;
  G4bool cutAsFinalRange;
  G4bool applyCuts;
  G4bool fluo;
  G4bool beardenFluoDir;
  G4bool auger;
  G4bool augerCascade;
  G4bool pixe;
  G4bool deexIgnoreCut;
  G4bool lateralDisplacement;
  G4bool muhadLateralDisplacement;
  G4bool latDisplacementBeyondSafety;
  G4bool useAngGeneratorForIonisation;
  G4bool useMottCorrection;
  G4bool integral;
  G4bool birks;

  G4double minSubRange;
  G4double minKinEnergy;
  G4double maxKinEnergy;
  G4double maxKinEnergyCSDA;
  G4double lowestElectronEnergy;
  G4double lowestMuHadEnergy;
  G4double linLossLimit;
  G4double bremsTh;
  G4double lambdaFactor;
  G4double factorForAngleLimit;
  G4double thetaLimit;
  G4double rangeFactor;
  G4double rangeFactorMuHad;
  G4double geomFactor;
  G4double skin;
  G4double dRoverRange;
  G4double finalRange;
  G4double dRoverRangeMuHad;
  G4double finalRangeMuHad;
  G4double factorScreen;

  G4int nbins;
  G4int nbinsPerDecade;
  G4int verbose;
  G4int workerVerbose;

  G4MscStepLimitType mscStepLimit;
  G4MscStepLimitType mscStepLimitMuHad;
  G4NuclearFormfactorType nucFormfactor;

  G4String namePIXE;
  G4String nameElectronPIXE;

  std::vector<G4String>  m_particlesPAI;
  std::vector<G4String>  m_regnamesPAI;
  std::vector<G4String>  m_typesPAI;

  std::vector<G4String>  m_regnamesME;

  std::vector<G4String>  m_regnamesDNA;
  std::vector<G4String>  m_typesDNA;

  std::vector<G4String>  m_regnamesMsc;
  std::vector<G4String>  m_typesMsc;

  std::vector<G4String>  m_regnamesSubCut;
  std::vector<G4bool>    m_subCuts;

  std::vector<G4String>  m_regnamesDeex;
  std::vector<G4bool>    m_fluo;
  std::vector<G4bool>    m_auger;
  std::vector<G4bool>    m_pixe;

  std::vector<G4String>  m_procBiasedXS;
  std::vector<G4double>  m_factBiasedXS;
  std::vector<G4bool>    m_weightBiasedXS;

  std::vector<G4String>  m_procForced;
  std::vector<G4String>  m_regnamesForced;
  std::vector<G4double>  m_lengthForced;
  std::vector<G4bool>    m_weightForced;

  std::vector<G4String>  m_procBiasedSec;
  std::vector<G4String>  m_regnamesBiasedSec;
  std::vector<G4double>  m_factBiasedSec;
  std::vector<G4double>  m_elimBiasedSec;

#ifdef G4MULTITHREADED
  static G4Mutex emParametersMutex;
#endif
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

