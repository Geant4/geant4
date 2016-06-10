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
// $Id: G4H2Messenger.hh 66310 2012-12-17 11:56:35Z ihrivnac $

// The messenger class for H2 management.
// It implements commands in /analysis/h2 directory.
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4H2Messenger_h
#define G4H2Messenger_h 1

#include "G4UImessenger.hh"
#include "G4AnalysisMessengerHelper.hh"
#include "globals.hh"

#include <memory>

class G4VAnalysisManager;
class G4UIdirectory;
class G4UIcommand;

class G4H2Messenger : public G4UImessenger
{
  public:
    explicit G4H2Messenger(G4VAnalysisManager* manager);
    virtual ~G4H2Messenger();
   
    // methods
    virtual void SetNewValue(G4UIcommand* command, G4String value) final;
    
  private:
    void CreateH2Cmd();
    void SetH2Cmd();
 
    G4VAnalysisManager*  fManager; ///< Associated class
    std::unique_ptr<G4AnalysisMessengerHelper>  fHelper;     
    std::unique_ptr<G4UIdirectory>  fDirectory;
    
    std::unique_ptr<G4UIcommand>  fCreateH2Cmd;
    std::unique_ptr<G4UIcommand>  fSetH2Cmd;
    std::unique_ptr<G4UIcommand>  fSetH2XCmd;
    std::unique_ptr<G4UIcommand>  fSetH2YCmd;
    std::unique_ptr<G4UIcommand>  fSetH2TitleCmd;   
    std::unique_ptr<G4UIcommand>  fSetH2XAxisCmd;   
    std::unique_ptr<G4UIcommand>  fSetH2YAxisCmd;   
    std::unique_ptr<G4UIcommand>  fSetH2ZAxisCmd;  

    G4int fXId;
    G4AnalysisMessengerHelper::BinData  fXData;
};
  
#endif
