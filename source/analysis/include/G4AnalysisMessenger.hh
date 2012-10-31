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
// $Id$

// Author: Ivana Hrivnacova, 04/07/2012  (ivana@ipno.in2p3.fr)
//
// This messenger class is a generalization of the HistoMessenger class,
// originally developed for the extended/electromagnetic examples
// by Michel Maire (michel.maire@lapp.in2p3.fr)

#ifndef G4AnalysisMessenger_h
#define G4AnalysisMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4VAnalysisManager;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;

class G4AnalysisMessenger : public G4UImessenger
{
  public:
    G4AnalysisMessenger(G4VAnalysisManager* manager);
    virtual ~G4AnalysisMessenger();
   
    // methods
    virtual void SetNewValue(G4UIcommand* command, G4String value);
    
  private:
    void CreateH1Cmd();
    void SetH1Cmd();
    void SetH1TitleCmd();
    void SetH1XAxisCmd();
    void SetH1YAxisCmd();
    void SetH1ActivationCmd();
    void CreateH2Cmd();
    void SetH2Cmd();
    void SetH2TitleCmd();
    void SetH2XAxisCmd();
    void SetH2YAxisCmd();
    void SetH2ZAxisCmd();
    void SetH2ActivationCmd();
 
    G4VAnalysisManager*  fManager; ///< Associated class
    
    G4UIdirectory*         fAnalysisDir;   
    G4UIcmdWithAString*    fSetFileNameCmd;
    G4UIcmdWithAString*    fSetHistoDirNameCmd;
    G4UIcmdWithAString*    fSetNtupleDirNameCmd;
    G4UIcmdWithABool*      fSetActivationCmd;   
    G4UIcmdWithAnInteger*  fVerboseCmd;   
    G4UIdirectory*         fH1Dir;   
    G4UIcommand*           fCreateH1Cmd;
    G4UIcommand*           fSetH1Cmd;
    G4UIcmdWithAnInteger*  fSetH1AsciiCmd;   
    G4UIcommand*           fSetH1TitleCmd;   
    G4UIcommand*           fSetH1XAxisCmd;   
    G4UIcommand*           fSetH1YAxisCmd;   
    G4UIcommand*           fSetH1ActivationCmd;   
    G4UIcmdWithABool*      fSetH1ActivationAllCmd;   
    G4UIdirectory*         fH2Dir;   
    G4UIcommand*           fCreateH2Cmd;
    G4UIcommand*           fSetH2Cmd;
    G4UIcmdWithAnInteger*  fSetH2AsciiCmd;   
    G4UIcommand*           fSetH2TitleCmd;   
    G4UIcommand*           fSetH2XAxisCmd;   
    G4UIcommand*           fSetH2YAxisCmd;   
    G4UIcommand*           fSetH2ZAxisCmd;   
    G4UIcommand*           fSetH2ActivationCmd;   
    G4UIcmdWithABool*      fSetH2ActivationAllCmd;   
};
  
#endif

