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
//
// $Id: Tst33AppStarterMessenger.hh,v 1.9 2008-04-21 09:00:03 ahoward Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class Tst33AppStarterMessenger
//
// Class description:
//
// Implementing the comands to message the application.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33AppStarterMessenger_hh
#define Tst33AppStarterMessenger_hh Tst33AppStarterMessenger_hh

#include <map>
#include "G4UImessenger.hh"

class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class Tst33AppStarter;
class Tst33VApplication;

typedef std::map<G4UIcmdWithAnInteger *, G4String> Tst33SimComands;

class Tst33AppStarterMessenger : public G4UImessenger{
public:
  explicit Tst33AppStarterMessenger(Tst33AppStarter &);
  virtual ~Tst33AppStarterMessenger();
  virtual void SetNewValue(G4UIcommand* pCmd,G4String szValue);

private:
  Tst33AppStarterMessenger(const Tst33AppStarterMessenger &);
  Tst33AppStarterMessenger &operator=(const Tst33AppStarterMessenger &);
  Tst33AppStarter &fAppStarter;
  G4UIcommand *fUseCoupledCmd;
  G4UIcommand *fMassGeoCmd;
  G4UIcommand *fParallelGeoCmd;
  G4UIcommand *fScoringCmd;
  G4UIcommand *fImpCmd;
  G4UIcmdWithAString *fWWCmd;
  G4UIcmdWithAnInteger *fWWRCmd;
  G4UIcommand *fClearSmaplingCmd;
  G4UIcommand *fConfigureSamplingCmd;
  G4UIcommand *fVisAppComand;
  G4UIcmdWithAnInteger *fTimedAppComand;
  G4UIcommand *fPostRunCmd;
  G4UIcmdWithAnInteger *fRunCmd;
  G4UIcommand *fWeightChangerCmd;
};

#endif
