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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRRunActionMessenger.hh
//   Header file of a messenger class that handles histograms and
//   n-tuple in GRRunAction.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#ifndef GRRunActionMessenger_H
#define GRRunActionMessenger_H 1

#include "G4UImessenger.hh"
#include "globals.hh"

class GRRunAction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;

class GRRunActionMessenger: public G4UImessenger
{
  public:
    GRRunActionMessenger(GRRunAction*);
    virtual ~GRRunActionMessenger();
    virtual void SetNewValue(G4UIcommand*,G4String);
    virtual G4String GetCurrentValue(G4UIcommand*);

  private:
    void Define1D();
    void Define1P();

  private:
    GRRunAction* pRA;

    G4UIdirectory* 		anaDir;
    G4UIcmdWithAString*         fileCmd;
    G4UIcmdWithAnInteger*	verboseCmd;
    G4UIcmdWithoutParameter*    listCmd;
    G4UIcmdWithAnInteger*	openCmd;
    G4UIcmdWithAnInteger*       plotCmd;
    G4UIcmdWithABool*           carryCmd;
    G4UIcmdWithoutParameter*    flushCmd;
    G4UIcmdWithoutParameter*    resetCmd;
    G4UIcommand*                idOffsetCmd;

    G4UIdirectory*              oneDDir;
    G4UIcommand*                create1DCmd;
    G4UIcommand*                create1DPrimPCmd;
    G4UIcommand*                create1DPlotPCmd;
    G4UIcmdWithoutParameter*    set1DCmd;
    G4UIcommand*                config1DCmd;
    G4UIcommand*                title1DCmd;
    G4UIcmdWithABool*           set1DYaxisLogCmd;

    G4UIdirectory*              onePDir;
    G4UIcommand*                create1PCmd;
    G4UIcommand*                set1PCmd;
    G4UIcommand*                title1PCmd;

    G4UIdirectory*              ntupleDir;
    G4UIcommand*                addColumnCmd;

  private:
    G4int currentID;

  private:
    inline G4bool CheckID(G4UIcommand* cmd)
    {
      if(currentID<0)
      {
        G4ExceptionDescription ed;
        ed << "There is no currently opened histogram. Create or open one.";
        cmd->CommandFailed(ed);
        return false;
      }
      return true;
    }
    inline G4bool CheckOpenID(G4UIcommand* /*cmd*/)
    {
      if(currentID>=0)
      {
        G4cout << "Previously opened histogram <" << currentID << "> is closed." << G4endl;
        currentID = -1;
      }
      return true;
    }
};

#endif

