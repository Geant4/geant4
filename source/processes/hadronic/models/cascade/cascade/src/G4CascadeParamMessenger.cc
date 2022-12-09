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
// Define simple UI commands as alternative to environment variables
//
// 20130304  M. Kelsey -- Add flag to collect and display cascade structure
// 20130621  M. Kelsey -- Add flag for CHECK_ECONS, replacing #ifdef's; add
//		flag to use three-body momentum parametrizations
// 20130703  M. Kelsey -- Add flag for USE_PHASESPACE
// 20141030  M. Kelsey -- Add flag to enable direct pi-N absorption
// 20141211  M. Kelsey -- Change PIN_ABSORPTION flag to double, for energy cut
// 20200110  M. Kelsey -- Reset cmdDir to 0 before .../cascade/ directory.

#include "G4CascadeParamMessenger.hh"
#include "G4CascadeParameters.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandTree.hh"
#include "G4UIdirectory.hh"
#include "G4UImanager.hh"


// Constructor and destructor

G4CascadeParamMessenger::G4CascadeParamMessenger(G4CascadeParameters* params)
  : G4UImessenger(), theParams(params), cmdDir(0), localCmdDir(false) {
  // NOTE: Put under same top-level tree as EM
  CreateDirectory("/process/had/","Hadronic processes"); cmdDir=0;
  CreateDirectory("/process/had/cascade/","Bertini-esque cascade parameters");

  verboseCmd = CreateCommand<G4UIcmdWithAnInteger>("verbose",
			"Enable information messages");
  balanceCmd = CreateCommand<G4UIcmdWithABool>("checkBalance",
			"Enable internal conservation checking");
  reportCmd = CreateCommand<G4UIcmdWithoutParameter>("report",
			"Dump all non-default parameter settings");
  usePreCoCmd = CreateCommand<G4UIcmdWithABool>("usePreCompound",
			"Use PreCompoundModel for nuclear de-excitation");
  doCoalCmd = CreateCommand<G4UIcmdWithABool>("doCoalescence",
			"Apply final-state nucleon clustering");
  piNAbsCmd = CreateCommand<G4UIcmdWithADouble>("piNAbsorption",
			"Probability for pion absorption on single nucleon");
  historyCmd = CreateCommand<G4UIcmdWithABool>("showHistory",
		        "Collect and report full structure of cascade");
  use3BodyCmd = CreateCommand<G4UIcmdWithABool>("use3BodyMom",
			"Use three-body momentum parametrizations");
  usePSCmd = CreateCommand<G4UIcmdWithABool>("usePhaseSpace",
			"Use Kopylov N-body momentum generator");
  randomFileCmd = CreateCommand<G4UIcmdWithAString>("randomFile",
			"Save random-engine to file at each interaction");
  nucUseBestCmd = CreateCommand<G4UIcmdWithABool>("useBestNuclearModel",
			"Use all physical-units for nuclear structure");
  nucRad2parCmd = CreateCommand<G4UIcmdWithADouble>("useTwoParamNuclearRadius",
			"Use R = C1*cbrt(A) + C2/cbrt(A)");
  nucRadScaleCmd = CreateCommand<G4UIcmdWithADouble>("nuclearRadiusScale",
			"Set length scale for nuclear model");
  nucRadSmallCmd = CreateCommand<G4UIcmdWithADouble>("smallNucleusRadius",
			"Set radius of A<4 nuclei");
  nucRadAlphaCmd = CreateCommand<G4UIcmdWithADouble>("alphaRadiusScale",
			"Fraction of small-radius for He-4");
  nucRadTrailingCmd = CreateCommand<G4UIcmdWithADouble>("shadowningRadius",
			"Effective nucleon radius for trailing effect");
  nucFermiScaleCmd = CreateCommand<G4UIcmdWithADouble>("fermiScale",
			"Scale factor for fermi momentum");
  nucXsecScaleCmd = CreateCommand<G4UIcmdWithADouble>("crossSectionScale",
			"Scale fator for total cross-sections");
  nucGammaQDCmd = CreateCommand<G4UIcmdWithADouble>("gammaQuasiDeutScale",
			"Scale factor for gamma-quasideutron cross-sections");
  coalDPmax2Cmd = CreateCommand<G4UIcmdWithADouble>("cluster2DPmax",
			"Maximum momentum for p-n clusters");
  coalDPmax3Cmd = CreateCommand<G4UIcmdWithADouble>("cluster3DPmax",
			"Maximum momentum for ppn/pnn clusters");
  coalDPmax4Cmd = CreateCommand<G4UIcmdWithADouble>("cluster4DPmax",
			"Maximum momentum for alpha clusters");
}

G4CascadeParamMessenger::~G4CascadeParamMessenger() {
  delete verboseCmd;
  delete balanceCmd;
  delete reportCmd;
  delete usePreCoCmd;
  delete doCoalCmd;
  delete piNAbsCmd;
  delete historyCmd;
  delete use3BodyCmd;
  delete usePSCmd;
  delete randomFileCmd;
  delete nucUseBestCmd;
  delete nucRad2parCmd;
  delete nucRadScaleCmd;
  delete nucRadSmallCmd;
  delete nucRadAlphaCmd;
  delete nucRadTrailingCmd;
  delete nucFermiScaleCmd;
  delete nucXsecScaleCmd;
  delete nucGammaQDCmd;
  delete coalDPmax2Cmd;
  delete coalDPmax3Cmd;
  delete coalDPmax4Cmd;
  if (localCmdDir) delete cmdDir;
}


// Create or reuse existing UIdirectory path

void G4CascadeParamMessenger::CreateDirectory(const char* path,
					      const char* desc) {
  G4UImanager* UIman = G4UImanager::GetUIpointer();
  if (!UIman) return;

  // Directory path must be absolute, prepend "/" if ncessary
  G4String fullPath = path;
  if (fullPath[0] != '/') fullPath.insert(0, '/', 1);
  if (fullPath.back() != '/') fullPath.append('/', 1);

  // See if input path has already been registered
  G4UIcommand* foundPath = UIman->GetTree()->FindPath(fullPath);
  if (foundPath) cmdDir = dynamic_cast<G4UIdirectory*>(foundPath);

  if (!cmdDir) {                // Create local deletable directory
    localCmdDir = true;
    cmdDir = new G4UIdirectory(fullPath.c_str());
    cmdDir->SetGuidance(desc);
  }
}


// Use command argument (literal string) to set envvar maps in container

void G4CascadeParamMessenger::SetNewValue(G4UIcommand* cmd, G4String arg) {
  if (cmd == reportCmd) theParams->DumpConfig(G4cout);

  if (cmd == verboseCmd) 
    theParams->G4CASCADE_VERBOSE = strdup(arg.c_str());

  if (cmd == balanceCmd) 
    theParams->G4CASCADE_CHECK_ECONS = StoB(arg) ? strdup(arg.c_str()) : 0;

  if (cmd == usePreCoCmd) 
    theParams->G4CASCADE_USE_PRECOMPOUND = StoB(arg) ? strdup(arg.c_str()) : 0;

  if (cmd == doCoalCmd) 
    theParams->G4CASCADE_DO_COALESCENCE = StoB(arg) ? strdup(arg.c_str()) : 0;

  if (cmd == piNAbsCmd) 
    theParams->G4CASCADE_PIN_ABSORPTION = strdup(arg.c_str());

  if (cmd == historyCmd) 
    theParams->G4CASCADE_SHOW_HISTORY = StoB(arg) ? strdup(arg.c_str()) : 0;

  if (cmd == use3BodyCmd)
    theParams->G4CASCADE_USE_3BODYMOM = StoB(arg) ? strdup(arg.c_str()) : 0;

  if (cmd == usePSCmd)
    theParams->G4CASCADE_USE_PHASESPACE = StoB(arg) ? strdup(arg.c_str()) : 0;

  if (cmd == randomFileCmd)
    theParams->G4CASCADE_RANDOM_FILE = arg.empty() ? 0 : strdup(arg.c_str());

  if (cmd == nucUseBestCmd)
    theParams->G4NUCMODEL_USE_BEST = StoB(arg) ? strdup(arg.c_str()) : 0;

//  if (cmd == nucRad2parCmd)
//    theParams->G4NUCMODEL_RAD_2PAR = strdup(arg.c_str());
  if (cmd == nucRad2parCmd)
    theParams->G4NUCMODEL_RAD_2PAR = StoB(arg) ? strdup(arg.c_str()) : 0;

  if (cmd == nucRadScaleCmd)
    theParams->G4NUCMODEL_RAD_SCALE = strdup(arg.c_str());

  if (cmd == nucRadSmallCmd)
    theParams->G4NUCMODEL_RAD_SMALL = strdup(arg.c_str());

  if (cmd == nucRadAlphaCmd)
    theParams->G4NUCMODEL_RAD_ALPHA = strdup(arg.c_str());

  if (cmd == nucRadTrailingCmd)
    theParams->G4NUCMODEL_RAD_TRAILING = strdup(arg.c_str());

  if (cmd == nucFermiScaleCmd)
    theParams->G4NUCMODEL_FERMI_SCALE = strdup(arg.c_str());

  if (cmd == nucXsecScaleCmd)
    theParams->G4NUCMODEL_XSEC_SCALE = strdup(arg.c_str());

  if (cmd == nucGammaQDCmd)
    theParams->G4NUCMODEL_GAMMAQD = strdup(arg.c_str());

  if (cmd == coalDPmax2Cmd)
    theParams->DPMAX_2CLUSTER = strdup(arg.c_str());

  if (cmd == coalDPmax3Cmd)
    theParams->DPMAX_3CLUSTER = strdup(arg.c_str());

  if (cmd == coalDPmax4Cmd)
    theParams->DPMAX_4CLUSTER = strdup(arg.c_str());

  theParams->Initialize();	// Update numerical values from settings
}
