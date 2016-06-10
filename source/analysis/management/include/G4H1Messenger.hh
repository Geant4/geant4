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
// $Id: G4H1Messenger.hh 66310 2012-12-17 11:56:35Z ihrivnac $

// The messenger class for H1 management.
// It implements commands in /analysis/h1 directory.
//
// This messenger class is a generalization of the HistoMessenger class,
// originally developed for the extended/electromagnetic examples
// by Michel Maire (michel.maire@lapp.in2p3.fr)
//
// Author: Ivana Hrivnacova, 24/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4H1Messenger_h
#define G4H1Messenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include <memory>

class G4VAnalysisManager;
class G4AnalysisMessengerHelper;
class G4UIdirectory;
class G4UIcommand;

class G4H1Messenger : public G4UImessenger
{
  public:
    explicit G4H1Messenger(G4VAnalysisManager* manager);
    virtual ~G4H1Messenger();
   
    // methods
    virtual void SetNewValue(G4UIcommand* command, G4String value) final;
    
  private:
    void CreateH1Cmd();
    void SetH1Cmd();
 
    G4VAnalysisManager*  fManager; ///< Associated class
    std::unique_ptr<G4AnalysisMessengerHelper>  fHelper; 
    std::unique_ptr<G4UIdirectory>  fDirectory;

    std::unique_ptr<G4UIcommand>  fCreateH1Cmd;
    std::unique_ptr<G4UIcommand>  fSetH1Cmd;
    std::unique_ptr<G4UIcommand>  fSetH1XCmd;
    std::unique_ptr<G4UIcommand>  fSetH1TitleCmd;   
    std::unique_ptr<G4UIcommand>  fSetH1XAxisCmd;   
    std::unique_ptr<G4UIcommand>  fSetH1YAxisCmd;   
};
  
#endif

