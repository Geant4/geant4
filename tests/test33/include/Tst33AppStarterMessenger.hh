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
// $Id: Tst33AppStarterMessenger.hh,v 1.2 2002-10-29 16:37:09 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class Tst33AppStarterMessenger
//
// Class description:
//
// ...

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33AppStarterMessenger_hh
#define Tst33AppStarterMessenger_hh Tst33AppStarterMessenger_hh

#include "g4std/map"
#include "G4UImessenger.hh"

class G4UIcmdWithAnInteger;
class Tst33AppStarter;
class Tst33VApplication;

typedef G4std::map<G4UIcmdWithAnInteger *, G4String> Tst33SimComands;

class Tst33AppStarterMessenger : public G4UImessenger{
public:
  explicit Tst33AppStarterMessenger(Tst33AppStarter &);
  virtual ~Tst33AppStarterMessenger();
  virtual void SetNewValue(G4UIcommand* pCmd,G4String szValue);

private:
  Tst33AppStarterMessenger(const Tst33AppStarterMessenger &);
  Tst33AppStarterMessenger &operator=(const Tst33AppStarterMessenger &);
  Tst33AppStarter &fAppStarter;
  G4UIcommand *fMassGeoCmd;
  G4UIcommand *fParallelGeoCmd;
  G4UIcommand *fScoringCmd;
  G4UIcommand *fImpCmd;
  G4UIcommand *fWWRCmd;
  G4UIcommand *fClearSmaplingCmd;
  G4UIcommand *fConfigureSamplingCmd;
  G4UIcommand *fVisAppComand;
  G4UIcmdWithAnInteger *fTimedAppComand;
  G4UIcommand *fPostRunCmd;
  G4UIcmdWithAnInteger *fRunCmd;
};

#endif
