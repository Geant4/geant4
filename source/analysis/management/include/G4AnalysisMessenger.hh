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
// $Id: G4AnalysisMessenger.hh 66310 2012-12-17 11:56:35Z ihrivnac $

// The messenger class for G4VAnalysisManager.
// It implements commands:
// - /analysis/setActivation
// - /analysis/verbose
//
// Author: Ivana Hrivnacova, 24/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4AnalysisMessenger_h
#define G4AnalysisMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4VAnalysisManager;
class G4HnManager;
class G4FileMessenger;
class G4H1Messenger;
class G4H2Messenger;
class G4H3Messenger;
class G4P1Messenger;
class G4P2Messenger;
class G4HnMessenger;

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;

class G4AnalysisMessenger : public G4UImessenger
{
  public:
    G4AnalysisMessenger(G4VAnalysisManager* manager);
    virtual ~G4AnalysisMessenger();
   
    // methods
    void SetH1HnManager(G4HnManager* h1HnManager);
    void SetH2HnManager(G4HnManager* h2HnManager);
    void SetH3HnManager(G4HnManager* h2HnManager);
    void SetP1HnManager(G4HnManager* h1HnManager);
    void SetP2HnManager(G4HnManager* h2HnManager);

    // methods
    virtual void SetNewValue(G4UIcommand* command, G4String value);

    // data members
    G4VAnalysisManager*    fManager; ///< Associated class
    G4FileMessenger*       fFileMessenger;
    G4H1Messenger*         fH1Messenger;
    G4H2Messenger*         fH2Messenger;
    G4H3Messenger*         fH3Messenger;
    G4P1Messenger*         fP1Messenger;
    G4P2Messenger*         fP2Messenger;
    G4HnMessenger*         fH1HnMessenger;
    G4HnMessenger*         fH2HnMessenger;
    G4HnMessenger*         fH3HnMessenger;
    G4HnMessenger*         fP1HnMessenger;
    G4HnMessenger*         fP2HnMessenger;
    
    G4UIdirectory*         fAnalysisDir;   
    G4UIcmdWithABool*      fSetActivationCmd;   
    G4UIcmdWithAnInteger*  fVerboseCmd;   
};
  
#endif

