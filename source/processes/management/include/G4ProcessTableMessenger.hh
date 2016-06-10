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
// $Id: G4ProcessTableMessenger.hh 71231 2013-06-12 13:06:28Z gcosmo $
//
//
//---------------------------------------------------------------
//
//  G4ProcessTableMessenger.hh
//
//  Class Description:
//    This is a messenger class to interface to exchange information
//    between ProcessTable and UI.
//-
//  /process/   Process Table control commands.
//   Commands : 
//     list *       :Dump process name registered.
//     verbose *    :Set Verbose Level for Process Table
//     activate *   :Activate process  
//     inactivate * :Inctivate process  
//     dump *       :Dump process information
//
//  History:
//    15 Aug. 1998, H. Kurashige  
//
//---------------------------------------------------------------

#ifndef G4ProcessTableMessenger_h
#define G4ProcessTableMessenger_h 1

class G4ProcessTable;

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger; 
class G4UIcmdWithAString;

#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4ProcessType.hh"

class G4ProcessTableMessenger: public G4UImessenger
{
  public:
    G4ProcessTableMessenger(G4ProcessTable* pTable);
    virtual ~G4ProcessTableMessenger();

  public: // with description
    virtual void SetNewValue(G4UIcommand * command,G4String newValues);
    virtual G4String GetCurrentValue(G4UIcommand * command);

  private:
    G4String GetProcessTypeName(G4ProcessType aType) const;  
    G4int    GetProcessType(const G4String& aTypeName) const;
    static G4ThreadLocal   G4int NumberOfProcessType;
    void     SetNumberOfProcessType();

  private:
    G4ProcessTableMessenger(const G4ProcessTableMessenger&)
      : G4UImessenger() {};
    G4ProcessTableMessenger(){};

  private:
    G4ProcessTable*  theProcessTable;

    G4UIdirectory *             thisDirectory;
    G4UIcmdWithAnInteger *      verboseCmd;
    G4UIcmdWithAString *        listCmd;
    G4UIcommand *               dumpCmd;
    G4UIcommand *               activateCmd;
    G4UIcommand *               inactivateCmd;
    G4UIcommand *               procVerboseCmd;

    G4String                    currentProcessTypeName;
    G4String                    currentProcessName;
    G4String                    currentParticleName; 
};

#endif




