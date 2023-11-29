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
// GEANT4 Class file
//
// File name:     G4EmParametersMessenger
//
// Author:        Vladimir Ivanchenko created from G4EnergyLossMessenger
//
// Creation date: 22-05-2013 
//
// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmParametersMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImanager.hh"
#include "G4MscStepLimitType.hh"
#include "G4NuclearFormfactorType.hh"
#include "G4EmParameters.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmParametersMessenger::G4EmParametersMessenger(G4EmParameters* ptr) 
  : theParameters(ptr)
{
  emDirectory = new G4UIdirectory("/process/em/", false);
  emDirectory->SetGuidance("General commands for EM processes.");
  eLossDirectory = new G4UIdirectory("/process/eLoss/", false);
  eLossDirectory->SetGuidance("Commands for energy loss processes.");
  mscDirectory = new G4UIdirectory("/process/msc/", false);
  mscDirectory->SetGuidance("Commands for EM scattering processes.");
  gconvDirectory = new G4UIdirectory("/process/gconv/", false);
  gconvDirectory->SetGuidance("Commands for EM gamma conversion BH5D model.");
  dnaDirectory = new G4UIdirectory("/process/dna/", false);
  dnaDirectory->SetGuidance("Commands for DNA processes.");

  flucCmd = new G4UIcmdWithABool("/process/eLoss/fluct",this);
  flucCmd->SetGuidance("Enable/disable energy loss fluctuations.");
  flucCmd->SetParameterName("choice",true);
  flucCmd->SetDefaultValue(true);
  flucCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  flucCmd->SetToBeBroadcasted(false);

  rangeCmd = new G4UIcmdWithABool("/process/eLoss/CSDARange",this);
  rangeCmd->SetGuidance("Enable/disable CSDA range calculation");
  rangeCmd->SetParameterName("range",true);
  rangeCmd->SetDefaultValue(false);
  rangeCmd->AvailableForStates(G4State_PreInit);
  rangeCmd->SetToBeBroadcasted(false);

  lpmCmd = new G4UIcmdWithABool("/process/eLoss/LPM",this);
  lpmCmd->SetGuidance("Enable/disable LPM effect calculation");
  lpmCmd->SetParameterName("lpm",true);
  lpmCmd->SetDefaultValue(true);
  lpmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  lpmCmd->SetToBeBroadcasted(false);

  rsCmd = new G4UIcmdWithABool("/process/eLoss/useCutAsFinalRange",this);
  rsCmd->SetGuidance("Enable/disable use of cut in range as a final range");
  rsCmd->SetParameterName("choice",true);
  rsCmd->SetDefaultValue(false);
  rsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  rsCmd->SetToBeBroadcasted(false);

  aplCmd = new G4UIcmdWithABool("/process/em/applyCuts",this);
  aplCmd->SetGuidance("Enable/disable applying cuts for gamma processes");
  aplCmd->SetParameterName("apl",true);
  aplCmd->SetDefaultValue(false);
  aplCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  aplCmd->SetToBeBroadcasted(false);

  intCmd = new G4UIcmdWithABool("/process/em/integral",this);
  intCmd->SetGuidance("Enable/disable integral method.");
  intCmd->SetParameterName("choice",true);
  intCmd->SetDefaultValue(true);
  intCmd->AvailableForStates(G4State_PreInit);
  intCmd->SetToBeBroadcasted(false);

  latCmd = new G4UIcmdWithABool("/process/msc/LateralDisplacement",this);
  latCmd->SetGuidance("Enable/disable sampling of lateral displacement");
  latCmd->SetParameterName("lat",true);
  latCmd->SetDefaultValue(true);
  latCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  latCmd->SetToBeBroadcasted(false);

  lat96Cmd = new G4UIcmdWithABool("/process/msc/LateralDisplacementAlg96",this);
  lat96Cmd->SetGuidance("Enable/disable sampling of lateral displacement");
  lat96Cmd->SetParameterName("lat96",true);
  lat96Cmd->SetDefaultValue(false);
  lat96Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  lat96Cmd->SetToBeBroadcasted(false);

  mulatCmd = new G4UIcmdWithABool("/process/msc/MuHadLateralDisplacement",this);
  mulatCmd->SetGuidance("Enable/disable sampling of lateral displacement for muons and hadrons");
  mulatCmd->SetParameterName("mulat",true);
  mulatCmd->SetDefaultValue(true);
  mulatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  mulatCmd->SetToBeBroadcasted(false);

  delCmd = new G4UIcmdWithABool("/process/eLoss/UseAngularGenerator",this);
  delCmd->SetGuidance("Enable usage of angular generator for ionisation");
  delCmd->SetParameterName("del",true);
  delCmd->SetDefaultValue(false);
  delCmd->AvailableForStates(G4State_PreInit);
  delCmd->SetToBeBroadcasted(false);

  mottCmd = new G4UIcmdWithABool("/process/msc/UseMottCorrection",this);
  mottCmd->SetGuidance("Enable usage of Mott corrections for e- elastic scattering");
  mottCmd->SetParameterName("mott",true);
  mottCmd->SetDefaultValue(false);
  mottCmd->AvailableForStates(G4State_PreInit);
  mottCmd->SetToBeBroadcasted(false);

  birksCmd = new G4UIcmdWithABool("/process/em/UseG4EmSaturation",this);
  birksCmd->SetGuidance("Enable usage of built-in Birks saturation");
  birksCmd->SetParameterName("birks",true);
  birksCmd->SetDefaultValue(false);
  birksCmd->AvailableForStates(G4State_PreInit,G4State_Init);
  birksCmd->SetToBeBroadcasted(false);

  sharkCmd = new G4UIcmdWithABool("/process/em/UseGeneralProcess",this);
  sharkCmd->SetGuidance("Enable gamma, e+- general process");
  sharkCmd->SetParameterName("gen",true);
  sharkCmd->SetDefaultValue(false);
  sharkCmd->AvailableForStates(G4State_PreInit);
  sharkCmd->SetToBeBroadcasted(false);

  poCmd = new G4UIcmdWithABool("/process/em/Polarisation",this);
  poCmd->SetGuidance("Enable polarisation");
  poCmd->AvailableForStates(G4State_PreInit);
  poCmd->SetToBeBroadcasted(false);

  sampleTCmd = new G4UIcmdWithABool("/process/em/enableSamplingTable",this);
  sampleTCmd->SetGuidance("Enable usage of sampling table for secondary generation");
  sampleTCmd->SetParameterName("sampleT",true);
  sampleTCmd->SetDefaultValue(false);
  sampleTCmd->AvailableForStates(G4State_PreInit);
  sampleTCmd->SetToBeBroadcasted(false);

  icru90Cmd = new G4UIcmdWithABool("/process/eLoss/UseICRU90",this);
  icru90Cmd->SetGuidance("Enable usage of ICRU90 stopping powers");
  icru90Cmd->SetParameterName("icru90",true);
  icru90Cmd->SetDefaultValue(false);
  icru90Cmd->AvailableForStates(G4State_PreInit);
  icru90Cmd->SetToBeBroadcasted(false);

  mudatCmd = new G4UIcmdWithABool("/process/em/MuDataFromFile",this);
  mudatCmd->SetGuidance("Enable usage of muon data from file");
  mudatCmd->SetParameterName("mudat",true);
  mudatCmd->SetDefaultValue(false);
  mudatCmd->AvailableForStates(G4State_PreInit);
  mudatCmd->SetToBeBroadcasted(false);

  peKCmd = new G4UIcmdWithABool("/process/em/PhotoeffectBelowKShell",this);
  peKCmd->SetGuidance("Enable sampling of photoeffect below K-shell");
  peKCmd->SetParameterName("peK",true);
  peKCmd->SetDefaultValue(true);
  peKCmd->AvailableForStates(G4State_PreInit);
  peKCmd->SetToBeBroadcasted(false);

  mscPCmd = new G4UIcmdWithABool("/process/msc/PositronCorrection",this);
  mscPCmd->SetGuidance("Enable msc positron correction");
  mscPCmd->SetParameterName("mscPC",true);
  mscPCmd->SetDefaultValue(true);
  mscPCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  mscPCmd->SetToBeBroadcasted(false);

  minEnCmd = new G4UIcmdWithADoubleAndUnit("/process/eLoss/minKinEnergy",this);
  minEnCmd->SetGuidance("Set the min kinetic energy for EM tables");
  minEnCmd->SetParameterName("emin",true);
  minEnCmd->SetUnitCategory("Energy");
  minEnCmd->AvailableForStates(G4State_PreInit);
  minEnCmd->SetToBeBroadcasted(false);

  maxEnCmd = new G4UIcmdWithADoubleAndUnit("/process/eLoss/maxKinEnergy",this);
  maxEnCmd->SetGuidance("Set the max kinetic energy for EM tables");
  maxEnCmd->SetParameterName("emax",true);
  maxEnCmd->SetUnitCategory("Energy");
  maxEnCmd->AvailableForStates(G4State_PreInit);
  maxEnCmd->SetToBeBroadcasted(false);

  cenCmd = new G4UIcmdWithADoubleAndUnit("/process/eLoss/maxKinEnergyCSDA",this);
  cenCmd->SetGuidance("Set the max kinetic energy for CSDA table");
  cenCmd->SetParameterName("emaxCSDA",true);
  cenCmd->SetUnitCategory("Energy");
  cenCmd->AvailableForStates(G4State_PreInit);
  cenCmd->SetToBeBroadcasted(false);

  max5DCmd = new G4UIcmdWithADoubleAndUnit("/process/em/max5DMuPairEnergy",this);
  max5DCmd->SetGuidance("Set the max kinetic energy for 5D muon pair production");
  max5DCmd->SetParameterName("emax5D",true);
  max5DCmd->SetUnitCategory("Energy");
  max5DCmd->AvailableForStates(G4State_PreInit);
  max5DCmd->SetToBeBroadcasted(false);

  lowEnCmd = new G4UIcmdWithADoubleAndUnit("/process/em/lowestElectronEnergy",this);
  lowEnCmd->SetGuidance("Set the lowest kinetic energy for e+-");
  lowEnCmd->SetParameterName("elow",true);
  lowEnCmd->SetUnitCategory("Energy");
  lowEnCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  lowEnCmd->SetToBeBroadcasted(false);

  lowhEnCmd = new G4UIcmdWithADoubleAndUnit("/process/em/lowestMuHadEnergy",this);
  lowhEnCmd->SetGuidance("Set the lowest kinetic energy for muons and hadrons");
  lowhEnCmd->SetParameterName("elowh",true);
  lowhEnCmd->SetUnitCategory("Energy");
  lowhEnCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  lowhEnCmd->SetToBeBroadcasted(false);

  lowEn3Cmd = new G4UIcmdWithADoubleAndUnit("/process/em/lowestTripletEnergy",this);
  lowEn3Cmd->SetGuidance("Set the lowest kinetic energy for triplet production");
  lowEn3Cmd->SetParameterName("elow3",true);
  lowEn3Cmd->SetUnitCategory("Energy");
  lowEn3Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  lowEn3Cmd->SetToBeBroadcasted(false);

  lllCmd = new G4UIcmdWithADouble("/process/eLoss/linLossLimit",this);
  lllCmd->SetGuidance("Set linearLossLimit parameter");
  lllCmd->SetParameterName("linlim",true);
  lllCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  lllCmd->SetToBeBroadcasted(false);

  brCmd = new G4UIcmdWithADoubleAndUnit("/process/eLoss/bremThreshold",this);
  brCmd->SetGuidance("Set e+- bremsstrahlung energy threshold");
  brCmd->SetParameterName("emaxBrem",true);
  brCmd->SetUnitCategory("Energy");
  brCmd->AvailableForStates(G4State_PreInit);
  brCmd->SetToBeBroadcasted(false);

  br1Cmd = new G4UIcmdWithADoubleAndUnit("/process/eLoss/bremMuHadThreshold",this);
  br1Cmd->SetGuidance("Set muon/hadron bremsstrahlung energy threshold");
  br1Cmd->SetParameterName("emaxMuHadBrem",true);
  br1Cmd->SetUnitCategory("Energy");
  br1Cmd->AvailableForStates(G4State_PreInit);
  br1Cmd->SetToBeBroadcasted(false);

  labCmd = new G4UIcmdWithADouble("/process/eLoss/LambdaFactor",this);
  labCmd->SetGuidance("Set lambdaFactor parameter for integral option");
  labCmd->SetParameterName("Fl",true);
  labCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  labCmd->SetToBeBroadcasted(false);

  mscfCmd = new G4UIcmdWithADouble("/process/msc/FactorForAngleLimit",this);
  mscfCmd->SetGuidance("Set factor for computation of a limit for -t (invariant transfer)");
  mscfCmd->SetParameterName("Fact",true);
  mscfCmd->SetRange("Fact>0");
  mscfCmd->SetDefaultValue(1.);
  mscfCmd->AvailableForStates(G4State_PreInit);
  mscfCmd->SetToBeBroadcasted(false);

  angCmd = new G4UIcmdWithADoubleAndUnit("/process/msc/ThetaLimit",this);
  angCmd->SetGuidance("Set the limit on the polar angle for msc and single scattering");
  angCmd->SetParameterName("theta",true);
  angCmd->SetUnitCategory("Angle");
  angCmd->AvailableForStates(G4State_PreInit);
  angCmd->SetToBeBroadcasted(false);

  msceCmd = new G4UIcmdWithADoubleAndUnit("/process/msc/EnergyLimit",this);
  msceCmd->SetGuidance("Set the upper energy limit for msc");
  msceCmd->SetParameterName("mscE",true);
  msceCmd->SetUnitCategory("Energy");
  msceCmd->AvailableForStates(G4State_PreInit);
  msceCmd->SetToBeBroadcasted(false);

  nielCmd = new G4UIcmdWithADoubleAndUnit("/process/em/MaxEnergyNIEL",this);
  nielCmd->SetGuidance("Set the upper energy limit for NIEL");
  nielCmd->SetParameterName("niel",true);
  nielCmd->SetUnitCategory("Energy");
  nielCmd->AvailableForStates(G4State_PreInit);
  nielCmd->SetToBeBroadcasted(false);

  frCmd = new G4UIcmdWithADouble("/process/msc/RangeFactor",this);
  frCmd->SetGuidance("Set RangeFactor for msc processes of e+-");
  frCmd->SetParameterName("Fr",true);
  frCmd->SetRange("Fr>0");
  frCmd->SetDefaultValue(0.04);
  frCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  frCmd->SetToBeBroadcasted(false);

  fr1Cmd = new G4UIcmdWithADouble("/process/msc/RangeFactorMuHad",this);
  fr1Cmd->SetGuidance("Set RangeFactor for msc processes of muons/hadrons");
  fr1Cmd->SetParameterName("Fr1",true);
  fr1Cmd->SetRange("Fr1>0");
  fr1Cmd->SetDefaultValue(0.2);
  fr1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fr1Cmd->SetToBeBroadcasted(false);

  fgCmd = new G4UIcmdWithADouble("/process/msc/GeomFactor",this);
  fgCmd->SetGuidance("Set GeomFactor parameter for msc processes");
  fgCmd->SetParameterName("Fg",true);
  fgCmd->SetRange("Fg>0");
  fgCmd->SetDefaultValue(2.5);
  fgCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fgCmd->SetToBeBroadcasted(false);
 
  skinCmd = new G4UIcmdWithADouble("/process/msc/Skin",this);
  skinCmd->SetGuidance("Set skin parameter for msc processes");
  skinCmd->SetParameterName("skin",true);
  skinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  skinCmd->SetToBeBroadcasted(false);

  screCmd = new G4UIcmdWithADouble("/process/msc/ScreeningFactor",this);
  screCmd->SetGuidance("Set screening factor");
  screCmd->SetParameterName("screen",true);
  screCmd->AvailableForStates(G4State_PreInit);
  screCmd->SetToBeBroadcasted(false);

  safCmd = new G4UIcmdWithADouble("/process/msc/SafetyFactor",this);
  safCmd->SetGuidance("Set safety factor");
  safCmd->SetParameterName("fsafe",true);
  safCmd->AvailableForStates(G4State_PreInit);
  safCmd->SetToBeBroadcasted(false);

  llimCmd = new G4UIcmdWithADoubleAndUnit("/process/msc/LambdaLimit",this);
  llimCmd->SetGuidance("Set the upper energy limit for NIEL");
  llimCmd->SetParameterName("ll",true);
  llimCmd->SetUnitCategory("Length");
  llimCmd->AvailableForStates(G4State_PreInit);
  llimCmd->SetToBeBroadcasted(false);

  amCmd = new G4UIcmdWithAnInteger("/process/em/binsPerDecade",this);
  amCmd->SetGuidance("Set number of bins per decade for EM tables");
  amCmd->SetParameterName("bins",true);
  amCmd->SetDefaultValue(7);
  amCmd->AvailableForStates(G4State_PreInit);
  amCmd->SetToBeBroadcasted(false);

  verCmd = new G4UIcmdWithAnInteger("/process/eLoss/verbose",this);
  verCmd->SetGuidance("Set verbose level for EM physics");
  verCmd->SetParameterName("verb",true);
  verCmd->SetDefaultValue(1);
  verCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  verCmd->SetToBeBroadcasted(false);

  ver1Cmd = new G4UIcmdWithAnInteger("/process/em/verbose",this);
  ver1Cmd->SetGuidance("Set verbose level for EM physics");
  ver1Cmd->SetParameterName("verb1",true);
  ver1Cmd->SetDefaultValue(1);
  ver1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  ver1Cmd->SetToBeBroadcasted(false);

  ver2Cmd = new G4UIcmdWithAnInteger("/process/em/workerVerbose",this);
  ver2Cmd->SetGuidance("Set worker verbose level for EM physics");
  ver2Cmd->SetParameterName("verb2",true);
  ver2Cmd->SetDefaultValue(0);
  ver2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  ver2Cmd->SetToBeBroadcasted(false);

  transWithMscCmd = new G4UIcmdWithAString("/process/em/transportationWithMsc",this);
  transWithMscCmd->SetGuidance("Enable/disable the G4TransportationWithMsc process");
  transWithMscCmd->SetParameterName("trans",true);
  transWithMscCmd->SetCandidates("Disabled Enabled MultipleSteps");
  transWithMscCmd->AvailableForStates(G4State_PreInit);
  transWithMscCmd->SetToBeBroadcasted(false);

  mscCmd = new G4UIcmdWithAString("/process/msc/StepLimit",this);
  mscCmd->SetGuidance("Set msc step limitation type");
  mscCmd->SetParameterName("StepLim",true);
  mscCmd->SetCandidates("Minimal UseSafety UseSafetyPlus UseDistanceToBoundary");
  mscCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  mscCmd->SetToBeBroadcasted(false);

  msc1Cmd = new G4UIcmdWithAString("/process/msc/StepLimitMuHad",this);
  msc1Cmd->SetGuidance("Set msc step limitation type for muons/hadrons");
  msc1Cmd->SetParameterName("StepLim1",true);
  msc1Cmd->SetCandidates("Minimal UseSafety UseSafetyPlus UseDistanceToBoundary");
  msc1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  msc1Cmd->SetToBeBroadcasted(false);

  dumpCmd = new G4UIcommand("/process/em/printParameters",this);
  dumpCmd->SetGuidance("Print all EM parameters.");
  dumpCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  dumpCmd->SetToBeBroadcasted(false);

  nffCmd = new G4UIcmdWithAString("/process/em/setNuclearFormFactor",this);
  nffCmd->SetGuidance("Define type of nuclear form-factor");
  nffCmd->SetParameterName("NucFF",true);
  nffCmd->SetCandidates("None Exponential Gaussian Flat");
  nffCmd->AvailableForStates(G4State_PreInit);
  nffCmd->SetToBeBroadcasted(false);

  ssCmd = new G4UIcmdWithAString("/process/em/setSingleScattering",this);
  ssCmd->SetGuidance("Define type of e+- single scattering model");
  ssCmd->SetParameterName("SS",true);
  ssCmd->SetCandidates("WVI Mott DPWA");
  ssCmd->AvailableForStates(G4State_PreInit);
  ssCmd->SetToBeBroadcasted(false);

  fluc1Cmd = new G4UIcmdWithAString("/process/eloss/setFluctModel",this);
  fluc1Cmd->SetGuidance("Define type of energy loss fluctuation model");
  fluc1Cmd->SetParameterName("Fluc1",true);
  fluc1Cmd->SetCandidates("Dummy Universal Urban");
  fluc1Cmd->AvailableForStates(G4State_PreInit);
  fluc1Cmd->SetToBeBroadcasted(false);

  tripletCmd = new G4UIcmdWithAnInteger("/process/gconv/conversionType",this);
  tripletCmd->SetGuidance("gamma conversion triplet/nuclear generation type:");
  tripletCmd->SetGuidance("0 - (default) both triplet and nuclear");
  tripletCmd->SetGuidance("1 - force nuclear");
  tripletCmd->SetGuidance("2 - force triplet");
  tripletCmd->SetParameterName("type",false);
  tripletCmd->SetRange("type >= 0 && type <= 2");
  tripletCmd->SetDefaultValue(0);
  tripletCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  tripletCmd->SetToBeBroadcasted(false);

  onIsolatedCmd = new G4UIcmdWithABool("/process/gconv/onIsolated",this);
  onIsolatedCmd->SetGuidance("Conversion on isolated charged particles");
  onIsolatedCmd->SetGuidance("false (default) : atomic electron screening");
  onIsolatedCmd->SetGuidance("true : conversion on isolated particles.");
  onIsolatedCmd->SetParameterName("flag",false);
  onIsolatedCmd->SetDefaultValue(false);
  onIsolatedCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  onIsolatedCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmParametersMessenger::~G4EmParametersMessenger()
{
  delete gconvDirectory;
  delete eLossDirectory;
  delete mscDirectory;
  delete emDirectory;
  delete dnaDirectory;

  delete flucCmd;
  delete rangeCmd;
  delete lpmCmd;
  delete rsCmd;
  delete aplCmd;
  delete intCmd;
  delete latCmd;
  delete lat96Cmd;
  delete mulatCmd;
  delete delCmd;
  delete mottCmd;
  delete birksCmd;
  delete sharkCmd;
  delete onIsolatedCmd;
  delete sampleTCmd;
  delete poCmd;
  delete icru90Cmd;
  delete mudatCmd;
  delete peKCmd;
  delete mscPCmd;

  delete minEnCmd;
  delete maxEnCmd;
  delete max5DCmd;
  delete cenCmd;
  delete lowEnCmd;
  delete lowhEnCmd;
  delete lowEn3Cmd;
  delete lllCmd;
  delete brCmd;
  delete br1Cmd;
  delete labCmd;
  delete mscfCmd;
  delete angCmd;
  delete msceCmd;
  delete nielCmd;
  delete frCmd;
  delete fr1Cmd;
  delete fgCmd;
  delete skinCmd;
  delete safCmd;
  delete llimCmd;
  delete screCmd;

  delete amCmd;
  delete verCmd;
  delete ver1Cmd;
  delete ver2Cmd;
  delete transWithMscCmd;
  delete tripletCmd;

  delete mscCmd;
  delete msc1Cmd;
  delete nffCmd;
  delete ssCmd;
  delete fluc1Cmd;

  delete dumpCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmParametersMessenger::SetNewValue(G4UIcommand* command, 
                                          G4String newValue)
{
  G4bool physicsModified = false;
  if (command == flucCmd) {
    theParameters->SetLossFluctuations(flucCmd->GetNewBoolValue(newValue));
    physicsModified = true;
  } else if (command == rangeCmd) {
    theParameters->SetBuildCSDARange(rangeCmd->GetNewBoolValue(newValue));
  } else if (command == lpmCmd) {
    theParameters->SetLPM(lpmCmd->GetNewBoolValue(newValue));
    physicsModified = true;
  } else if (command == rsCmd) {
    theParameters->SetUseCutAsFinalRange(rsCmd->GetNewBoolValue(newValue));
    physicsModified = true;
  } else if (command == aplCmd) {
    theParameters->SetApplyCuts(aplCmd->GetNewBoolValue(newValue));
    physicsModified = true;
  } else if (command == intCmd) {
    theParameters->SetIntegral(intCmd->GetNewBoolValue(newValue));
  } else if (command == latCmd) {
    theParameters->SetLateralDisplacement(latCmd->GetNewBoolValue(newValue));
    physicsModified = true;
  } else if (command == lat96Cmd) {
    theParameters->SetLateralDisplacementAlg96(lat96Cmd->GetNewBoolValue(newValue));
    physicsModified = true;
  } else if (command == mulatCmd) {
    theParameters->SetMuHadLateralDisplacement(mulatCmd->GetNewBoolValue(newValue));
    physicsModified = true;
  } else if (command == delCmd) {
    theParameters->ActivateAngularGeneratorForIonisation(delCmd->GetNewBoolValue(newValue));
  } else if (command == mottCmd) {
    theParameters->SetUseMottCorrection(mottCmd->GetNewBoolValue(newValue));
  } else if (command == birksCmd) {
    theParameters->SetBirksActive(birksCmd->GetNewBoolValue(newValue));
  } else if (command == icru90Cmd) {
    theParameters->SetUseICRU90Data(icru90Cmd->GetNewBoolValue(newValue));
  } else if (command == sharkCmd) {
    theParameters->SetGeneralProcessActive(sharkCmd->GetNewBoolValue(newValue));
  } else if (command == poCmd) {
    theParameters->SetEnablePolarisation(poCmd->GetNewBoolValue(newValue));
  } else if (command == sampleTCmd) {
    theParameters->SetEnableSamplingTable(sampleTCmd->GetNewBoolValue(newValue));
  } else if (command == mudatCmd) {
    theParameters->SetRetrieveMuDataFromFile(mudatCmd->GetNewBoolValue(newValue));
  } else if (command == peKCmd) {
    theParameters->SetPhotoeffectBelowKShell(peKCmd->GetNewBoolValue(newValue));
  } else if (command == mscPCmd) {
    theParameters->SetMscPositronCorrection(mscPCmd->GetNewBoolValue(newValue));

  } else if (command == minEnCmd) {
    theParameters->SetMinEnergy(minEnCmd->GetNewDoubleValue(newValue));
  } else if (command == maxEnCmd) { 
    theParameters->SetMaxEnergy(maxEnCmd->GetNewDoubleValue(newValue));
  } else if (command == max5DCmd) { 
    theParameters->SetMaxEnergyFor5DMuPair(max5DCmd->GetNewDoubleValue(newValue));
  } else if (command == cenCmd) { 
    theParameters->SetMaxEnergyForCSDARange(cenCmd->GetNewDoubleValue(newValue));
    physicsModified = true;
  } else if (command == lowEnCmd) { 
    theParameters->SetLowestElectronEnergy(lowEnCmd->GetNewDoubleValue(newValue));
    physicsModified = true;
  } else if (command == lowEn3Cmd) { 
    theParameters->SetLowestTripletEnergy(lowEn3Cmd->GetNewDoubleValue(newValue));
    physicsModified = true;
  } else if (command == lowhEnCmd) { 
    theParameters->SetLowestMuHadEnergy(lowhEnCmd->GetNewDoubleValue(newValue));
    physicsModified = true;
  } else if (command == lllCmd) { 
    theParameters->SetLinearLossLimit(lllCmd->GetNewDoubleValue(newValue));
    physicsModified = true;
  } else if (command == brCmd) { 
    theParameters->SetBremsstrahlungTh(brCmd->GetNewDoubleValue(newValue));
    physicsModified = true;
  } else if (command == br1Cmd) { 
    theParameters->SetMuHadBremsstrahlungTh(br1Cmd->GetNewDoubleValue(newValue));
    physicsModified = true;
  } else if (command == labCmd) {
    theParameters->SetLambdaFactor(labCmd->GetNewDoubleValue(newValue));
    physicsModified = true;
  } else if (command == mscfCmd) {
    theParameters->SetFactorForAngleLimit(mscfCmd->GetNewDoubleValue(newValue));
  } else if (command == angCmd) { 
    theParameters->SetMscThetaLimit(angCmd->GetNewDoubleValue(newValue));
  } else if (command == msceCmd) { 
    theParameters->SetMscEnergyLimit(msceCmd->GetNewDoubleValue(newValue));
  } else if (command == nielCmd) { 
    theParameters->SetMaxNIELEnergy(nielCmd->GetNewDoubleValue(newValue));
  } else if (command == frCmd) {
    theParameters->SetMscRangeFactor(frCmd->GetNewDoubleValue(newValue));
    physicsModified = true;
  } else if (command == fr1Cmd) {
    theParameters->SetMscMuHadRangeFactor(fr1Cmd->GetNewDoubleValue(newValue));
    physicsModified = true;
  } else if (command == fgCmd) {
    theParameters->SetMscGeomFactor(fgCmd->GetNewDoubleValue(newValue));
    physicsModified = true;
  } else if (command == skinCmd) { 
    theParameters->SetMscSkin(skinCmd->GetNewDoubleValue(newValue));
    physicsModified = true;
  } else if (command == safCmd) { 
    theParameters->SetMscSafetyFactor(safCmd->GetNewDoubleValue(newValue));
  } else if (command == llimCmd) { 
    theParameters->SetMscLambdaLimit(llimCmd->GetNewDoubleValue(newValue));
  } else if (command == screCmd) { 
    theParameters->SetScreeningFactor(screCmd->GetNewDoubleValue(newValue));
  } else if (command == amCmd) {
    theParameters->SetNumberOfBinsPerDecade(amCmd->GetNewIntValue(newValue));
  } else if (command == verCmd) {
    theParameters->SetVerbose(verCmd->GetNewIntValue(newValue));
  } else if (command == ver1Cmd) {
    theParameters->SetVerbose(ver1Cmd->GetNewIntValue(newValue));
  } else if (command == ver2Cmd) {
    theParameters->SetWorkerVerbose(ver2Cmd->GetNewIntValue(newValue));
  } else if (command == dumpCmd) {
    theParameters->SetIsPrintedFlag(false);
    theParameters->Dump();
  } else if (command == transWithMscCmd) {
    G4TransportationWithMscType type = G4TransportationWithMscType::fDisabled;
    if(newValue == "Disabled") {
      type = G4TransportationWithMscType::fDisabled;
    } else if(newValue == "Enabled") {
      type = G4TransportationWithMscType::fEnabled;
    } else if(newValue == "MultipleSteps") {
      type = G4TransportationWithMscType::fMultipleSteps;
    } else {
      G4ExceptionDescription ed;
      ed << " TransportationWithMsc type <" << newValue << "> unknown!";
      G4Exception("G4EmParametersMessenger", "em0044", JustWarning, ed);
    }
    theParameters->SetTransportationWithMsc(type);
  } else if (command == mscCmd || command == msc1Cmd) {
    G4MscStepLimitType msctype = fUseSafety;
    if(newValue == "Minimal") { 
      msctype = fMinimal;  
    } else if(newValue == "UseDistanceToBoundary") { 
      msctype = fUseDistanceToBoundary;
    } else if(newValue == "UseSafety") { 
      msctype = fUseSafety; 
    } else if(newValue == "UseSafetyPlus") {
      msctype = fUseSafetyPlus; 
    } else {
      G4ExceptionDescription ed;
      ed << " StepLimit type <" << newValue << "> unknown!"; 
      G4Exception("G4EmParametersMessenger", "em0044", JustWarning, ed);
      return;
    }
    if (command == mscCmd) {
      theParameters->SetMscStepLimitType(msctype);
    } else {
      theParameters->SetMscMuHadStepLimitType(msctype);
    }
    physicsModified = true;
  } else if (command == nffCmd) {
    G4NuclearFormfactorType x = fNoneNF;
    if(newValue == "Exponential") { x = fExponentialNF; }
    else if(newValue == "Gaussian") { x = fGaussianNF; }
    else if(newValue == "Flat") { x = fFlatNF; }
    else if(newValue != "None") { 
      G4ExceptionDescription ed;
      ed << " NuclearFormFactor type <" << newValue << "> unknown!"; 
      G4Exception("G4EmParametersMessenger", "em0044", JustWarning, ed);
      return; 
    }
    theParameters->SetNuclearFormfactorType(x);
  } else if (command == ssCmd) {
    G4eSingleScatteringType x = fWVI;
    if(newValue == "DPWA") { x = fDPWA; }
    else if(newValue == "Mott") { x = fMott; }
    else if(newValue != "WVI") { 
      G4ExceptionDescription ed;
      ed << " G4eSingleScatteringType type <" << newValue << "> unknown!"; 
      G4Exception("G4EmParametersMessenger", "em0044", JustWarning, ed);
      return; 
    }
    theParameters->SetSingleScatteringType(x);
  } else if (command == fluc1Cmd) {
    G4EmFluctuationType x = fUniversalFluctuation;
    if(newValue == "Dummy") { x = fDummyFluctuation; }
    else if(newValue == "Urban") { x = fUrbanFluctuation; }
    theParameters->SetFluctuationType(x);
  } else if ( command==tripletCmd ) {
    theParameters->SetConversionType(tripletCmd->GetNewIntValue(newValue));
  } else if ( command==onIsolatedCmd ) {
    theParameters->SetOnIsolated(onIsolatedCmd->GetNewBoolValue(newValue));
    physicsModified = true;
  }
  
  if(physicsModified) {
    G4UImanager::GetUIpointer()->ApplyCommand("/run/physicsModified");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
