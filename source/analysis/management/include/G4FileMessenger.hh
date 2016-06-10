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
// $Id: G4FileMessenger.hh 66310 2012-12-17 11:56:35Z ihrivnac $

// The messenger class for File management
//
// It implements commands:
// - /analysis/setFileName name
// - /analysis/setHistoDirName name
// - /analysis/setNtupleDirName name
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4FileMessenger_h
#define G4FileMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include <memory>

class G4VAnalysisManager;
class G4UIcmdWithAString;

class G4FileMessenger : public G4UImessenger
{
  public:
    explicit G4FileMessenger(G4VAnalysisManager* manager);
    virtual ~G4FileMessenger();
   
    // methods
    virtual void SetNewValue(G4UIcommand* command, G4String value) final;
 
    G4VAnalysisManager*  fManager; ///< Associated class
    
    std::unique_ptr<G4UIcmdWithAString>  fSetFileNameCmd;
    std::unique_ptr<G4UIcmdWithAString>  fSetHistoDirNameCmd;
    std::unique_ptr<G4UIcmdWithAString>  fSetNtupleDirNameCmd;
};
  
#endif

