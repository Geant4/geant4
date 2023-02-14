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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
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
#include "G4DNAModelSubType.hh"
#include "G4EmFluoDirectory.hh"
#include "G4EmSaturation.hh"
#include "G4ThreeVector.hh"
#include <vector>

enum G4eSingleScatteringType
{
  fWVI = 0,
  fMott,
  fDPWA
};

enum class G4TransportationWithMscType
{
  fDisabled = 0,
  fEnabled,
  fMultipleSteps,
};

enum G4EmFluctuationType 
{
  fDummyFluctuation = 0,
  fUniversalFluctuation,
  fUrbanFluctuation
};

class G4EmParametersMessenger;
class G4EmExtraParameters;
class G4EmLowEParameters;
class G4VAtomDeexcitation;
class G4VEnergyLossProcess;
class G4VEmProcess;
class G4StateManager;

class G4EmParameters
{
public:

  static G4EmParameters* Instance();

  ~G4EmParameters();

  void SetDefaults();

  // printing
  void StreamInfo(std::ostream& os) const;
  void Dump();
  friend std::ostream& operator<< (std::ostream& os, const G4EmParameters&);

  // boolean flags
  void SetLossFluctuations(G4bool val);
  G4bool LossFluctuation() const;

  void SetBuildCSDARange(G4bool val);
  G4bool BuildCSDARange() const;

  void SetLPM(G4bool val);
  G4bool LPM() const;

  void SetUseCutAsFinalRange(G4bool val);
  G4bool UseCutAsFinalRange() const;

  void SetApplyCuts(G4bool val);
  G4bool ApplyCuts() const;

  void SetFluo(G4bool val);
  G4bool Fluo() const;

  G4EmFluoDirectory FluoDirectory() const;

  void SetFluoDirectory(G4EmFluoDirectory);
  void SetBeardenFluoDir(G4bool val);
  void SetANSTOFluoDir(G4bool val);
  void SetXDB_EADLFluoDir(G4bool val);

  G4bool BeardenFluoDir();
  G4bool ANSTOFluoDir();

  void SetAuger(G4bool val);
  void SetAugerCascade(G4bool val) { SetAuger(val); };
  G4bool Auger() const;
  G4bool AugerCascade() const { return Auger(); }

  void SetPixe(G4bool val);
  G4bool Pixe() const;

  void SetDeexcitationIgnoreCut(G4bool val);
  G4bool DeexcitationIgnoreCut() const;

  void SetLateralDisplacement(G4bool val);
  G4bool LateralDisplacement() const;

  void SetLateralDisplacementAlg96(G4bool val);
  G4bool LateralDisplacementAlg96() const;

  void SetMuHadLateralDisplacement(G4bool val);
  G4bool MuHadLateralDisplacement() const;

  void ActivateAngularGeneratorForIonisation(G4bool val);
  G4bool UseAngularGeneratorForIonisation() const;

  void SetUseMottCorrection(G4bool val);
  G4bool UseMottCorrection() const;

  void SetIntegral(G4bool val);
  G4bool Integral() const;

  void SetBirksActive(G4bool val);
  G4bool BirksActive() const;

  void SetUseICRU90Data(G4bool val);
  G4bool UseICRU90Data() const;

  void SetFluctuationType(G4EmFluctuationType val);
  G4EmFluctuationType FluctuationType() const;

  void SetDNAFast(G4bool val);
  G4bool DNAFast() const;

  void SetDNAStationary(G4bool val);
  G4bool DNAStationary() const;

  void SetDNAElectronMsc(G4bool val);
  G4bool DNAElectronMsc() const;

  // if general interaction is enabled then 
  // force interaction options should be disabled
  void SetGeneralProcessActive(G4bool val);
  G4bool GeneralProcessActive() const;

  void SetEnableSamplingTable(G4bool val);
  G4bool EnableSamplingTable() const;

  void SetEnablePolarisation(G4bool val);
  G4bool EnablePolarisation() const;

  G4bool GetDirectionalSplitting() const;
  void SetDirectionalSplitting(G4bool v);

  G4bool QuantumEntanglement() const;
  void SetQuantumEntanglement(G4bool v);

  G4bool RetrieveMuDataFromFile() const;
  void SetRetrieveMuDataFromFile(G4bool v);

  G4bool PhotoeffectBelowKShell() const;
  void SetPhotoeffectBelowKShell(G4bool v);

  G4bool MscPositronCorrection() const;
  void SetMscPositronCorrection(G4bool v);

  // 5d
  void SetOnIsolated(G4bool val);
  G4bool OnIsolated() const;

  void ActivateDNA();
  void SetIsPrintedFlag(G4bool val);
  G4bool IsPrintLocked() const;

  // double parameters with values
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

  void SetLowestTripletEnergy(G4double val);
  G4double LowestTripletEnergy() const;

  void SetLinearLossLimit(G4double val);
  G4double LinearLossLimit() const;

  void SetBremsstrahlungTh(G4double val);
  G4double BremsstrahlungTh() const;
  void SetMuHadBremsstrahlungTh(G4double val);
  G4double MuHadBremsstrahlungTh() const;

  void SetLambdaFactor(G4double val);
  G4double LambdaFactor() const;

  void SetFactorForAngleLimit(G4double val);
  G4double FactorForAngleLimit() const;

  void SetMscThetaLimit(G4double val);
  G4double MscThetaLimit() const;

  void SetMscEnergyLimit(G4double val);
  G4double MscEnergyLimit() const;

  void SetMscRangeFactor(G4double val);
  G4double MscRangeFactor() const;

  void SetMscMuHadRangeFactor(G4double val);
  G4double MscMuHadRangeFactor() const;

  void SetMscGeomFactor(G4double val);
  G4double MscGeomFactor() const;

  void SetMscSafetyFactor(G4double val);
  G4double MscSafetyFactor() const;

  void SetMscLambdaLimit(G4double val);
  G4double MscLambdaLimit() const;

  void SetMscSkin(G4double val);
  G4double MscSkin() const;

  void SetScreeningFactor(G4double val);
  G4double ScreeningFactor() const;

  void SetMaxNIELEnergy(G4double val);
  G4double MaxNIELEnergy() const;

  void SetMaxEnergyFor5DMuPair(G4double val);
  G4double MaxEnergyFor5DMuPair() const;

  void SetStepFunction(G4double v1, G4double v2);
  void SetStepFunctionMuHad(G4double v1, G4double v2);
  void SetStepFunctionLightIons(G4double v1, G4double v2);
  void SetStepFunctionIons(G4double v1, G4double v2);
  void FillStepFunction(const G4ParticleDefinition*, G4VEnergyLossProcess*) const;

  void SetDirectionalSplittingRadius(G4double r);
  G4double GetDirectionalSplittingRadius();

  void SetDirectionalSplittingTarget(const G4ThreeVector& v);
  G4ThreeVector GetDirectionalSplittingTarget() const;

  // integer parameters 
  
  void SetNumberOfBinsPerDecade(G4int val);
  G4int NumberOfBinsPerDecade() const;
  G4int NumberOfBins() const;

  void SetVerbose(G4int val);
  G4int Verbose() const;

  void SetWorkerVerbose(G4int val);
  G4int WorkerVerbose() const;

  void SetTransportationWithMsc(G4TransportationWithMscType val);
  G4TransportationWithMscType TransportationWithMsc() const;

  void SetMscStepLimitType(G4MscStepLimitType val);
  G4MscStepLimitType MscStepLimitType() const;

  void SetMscMuHadStepLimitType(G4MscStepLimitType val);
  G4MscStepLimitType MscMuHadStepLimitType() const;

  void SetSingleScatteringType(G4eSingleScatteringType val); 
  G4eSingleScatteringType SingleScatteringType() const;

  void SetNuclearFormfactorType(G4NuclearFormfactorType val);
  G4NuclearFormfactorType NuclearFormfactorType() const;

  void SetDNAeSolvationSubType(G4DNAModelSubType val);
  G4DNAModelSubType DNAeSolvationSubType() const;

  //5d
  void  SetConversionType(G4int val);
  G4int GetConversionType() const;

  // string parameters 
  void SetPIXECrossSectionModel(const G4String&);
  const G4String& PIXECrossSectionModel();

  void SetPIXEElectronCrossSectionModel(const G4String&);
  const G4String& PIXEElectronCrossSectionModel();

  void SetLivermoreDataDir(const G4String&);
  const G4String& LivermoreDataDir();

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

  void AddPhysics(const G4String& region, const G4String& type);
  const std::vector<G4String>& RegionsPhysics() const;
  const std::vector<G4String>& TypesPhysics() const;

  void SetSubCutRegion(const G4String& region = "");

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

  // define external saturation class
  void SetEmSaturation(G4EmSaturation*);
  // create and access saturation class
  G4EmSaturation* GetEmSaturation();

  // initialisation methods
  void DefineRegParamForLoss(G4VEnergyLossProcess*) const;
  void DefineRegParamForEM(G4VEmProcess*) const;
  void DefineRegParamForDeex(G4VAtomDeexcitation*) const;

  G4EmParameters(G4EmParameters &) = delete;
  G4EmParameters & operator=(const G4EmParameters &right) = delete;  

private:

  G4EmParameters();

  void Initialise();

  G4bool IsLocked() const;

  void PrintWarning(G4ExceptionDescription& ed) const; 

  static G4EmParameters* theInstance;

  G4EmParametersMessenger* theMessenger;
  G4EmExtraParameters* fBParameters;
  G4EmLowEParameters* fCParameters;
  G4StateManager*  fStateManager;
  G4EmSaturation*  emSaturation;

  G4bool lossFluctuation;
  G4bool buildCSDARange;
  G4bool flagLPM;
  G4bool cutAsFinalRange;
  G4bool applyCuts;
  G4bool lateralDisplacement;
  G4bool lateralDisplacementAlg96;
  G4bool muhadLateralDisplacement;
  G4bool useAngGeneratorForIonisation;
  G4bool useMottCorrection;
  G4bool integral;
  G4bool birks;
  G4bool fICRU90;
  G4bool gener;
  G4bool fSamplingTable;
  G4bool fPolarisation;
  G4bool fMuDataFromFile;
  G4bool fPEKShell;
  G4bool fMscPosiCorr;
  G4bool onIsolated; // 5d model conversion on free ions
  G4bool fDNA;
  G4bool fIsPrinted;
  
  G4double minKinEnergy;
  G4double maxKinEnergy;
  G4double maxKinEnergyCSDA;
  G4double max5DEnergyForMuPair;
  G4double lowestElectronEnergy;
  G4double lowestMuHadEnergy;
  G4double lowestTripletEnergy;
  G4double linLossLimit;
  G4double bremsTh;
  G4double bremsMuHadTh;
  G4double lambdaFactor;
  G4double factorForAngleLimit;
  G4double thetaLimit;
  G4double energyLimit;
  G4double maxNIELEnergy;
  G4double rangeFactor;
  G4double rangeFactorMuHad;
  G4double geomFactor;
  G4double safetyFactor;
  G4double lambdaLimit;
  G4double skin;
  G4double factorScreen;

  G4int nbinsPerDecade;
  G4int verbose;
  G4int workerVerbose;
  G4int tripletConv;  // 5d model triplet generation type

  G4TransportationWithMscType fTransportationWithMsc;
  G4MscStepLimitType mscStepLimit;
  G4MscStepLimitType mscStepLimitMuHad;
  G4NuclearFormfactorType nucFormfactor;
  G4eSingleScatteringType fSStype;
  G4EmFluctuationType fFluct;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

