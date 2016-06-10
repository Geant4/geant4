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
// $Id: G4CascadeParamMessenger.hh 72016 2013-07-03 16:24:15Z mkelsey $
// Define simple UI commands as alternative to environment variables
//
// 20130304  M. Kelsey -- Add flag to collect and display cascade structure
// 20130621  M. Kelsey -- Add flag for CHECK_ECONS, replacing #ifdef's; add
//		flag to use three-body momentum parametrizations
// 20130703  M. Kelsey -- Add flag for USE_PHASESPACE
// 20141030  M. Kelsey -- Add flag to enable direct pi-N absorption
// 20141211  M. Kelsey -- Change PIN_ABSORPTION flag to double, for energy cut

#ifndef G4CascadeParamMessenger_hh
#define G4CascadeParamMessenger_hh

#include "G4UImessenger.hh"
#include "globals.hh"

class G4CascadeParameters;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcommand;
class G4UIdirectory;


class G4CascadeParamMessenger : public G4UImessenger {
public:
  G4CascadeParamMessenger(G4CascadeParameters* params);
  virtual ~G4CascadeParamMessenger();

  // Interface command needed by G4UImanager -- subclasses should call back!
  virtual void SetNewValue(G4UIcommand* command, G4String newValue);

protected:		// These commands are intended to be in a base class
  void CreateDirectory(const char* path, const char* desc);

  template <class T>
  T* CreateCommand(const G4String& cmd, const G4String& desc);

private:
  G4CascadeParameters*  theParams;
  G4UIdirectory* cmdDir;
  G4bool localCmdDir;		// Flag if created vs. found directory

  G4UIcmdWithAnInteger* verboseCmd;
  G4UIcmdWithoutParameter* reportCmd;
  G4UIcmdWithABool*	balanceCmd;
  G4UIcmdWithABool*     usePreCoCmd;
  G4UIcmdWithABool*     doCoalCmd;
  G4UIcmdWithADouble*   piNAbsCmd;
  G4UIcmdWithABool*     historyCmd;
  G4UIcmdWithABool*     use3BodyCmd;
  G4UIcmdWithABool*     usePSCmd;
  G4UIcmdWithAString*   randomFileCmd;
  G4UIcmdWithABool*     nucUseBestCmd;
  G4UIcmdWithADouble*   nucRad2parCmd;
  G4UIcmdWithADouble*   nucRadScaleCmd;
  G4UIcmdWithADouble*   nucRadSmallCmd;
  G4UIcmdWithADouble*   nucRadAlphaCmd;
  G4UIcmdWithADouble*   nucRadTrailingCmd;
  G4UIcmdWithADouble*   nucFermiScaleCmd;
  G4UIcmdWithADouble*   nucXsecScaleCmd;
  G4UIcmdWithADouble*   nucGammaQDCmd;
  G4UIcmdWithADouble*   coalDPmax2Cmd;
  G4UIcmdWithADouble*   coalDPmax3Cmd;
  G4UIcmdWithADouble*   coalDPmax4Cmd;
};

// Templated function implementation below
#include "G4CascadeParamMessenger.icc"

#endif /* G4CascadeParamMessenger_hh */
