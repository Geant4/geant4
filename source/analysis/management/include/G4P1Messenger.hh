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

// The messenger class for P1 management.
// It implements commands in /analysis/p1 directory.
//
// Author: Ivana Hrivnacova, 24/07/2014  (ivana@ipno.in2p3.fr)

#ifndef G4P1Messenger_h
#define G4P1Messenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4VAnalysisManager;
class G4UIdirectory;
class G4UIcommand;

class G4P1Messenger : public G4UImessenger
{
  public:
    G4P1Messenger(G4VAnalysisManager* manager);
    virtual ~G4P1Messenger();
   
    // methods
    virtual void SetNewValue(G4UIcommand* command, G4String value);
    
  private:
    void CreateP1Cmd();
    void SetP1Cmd();
    void SetP1TitleCmd();
    void SetP1XAxisCmd();
    void SetP1YAxisCmd();
 
    G4VAnalysisManager*  fManager; ///< Associated class
    
    G4UIdirectory*         fP1Dir;   
    G4UIcommand*           fCreateP1Cmd;
    G4UIcommand*           fSetP1Cmd;
    G4UIcommand*           fSetP1TitleCmd;   
    G4UIcommand*           fSetP1XAxisCmd;   
    G4UIcommand*           fSetP1YAxisCmd;   
};
  
#endif
