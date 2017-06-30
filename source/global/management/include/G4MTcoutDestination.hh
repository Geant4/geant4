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
// $Id: G4MTcoutDestination.hh 66241 2012-12-13 18:34:42Z gunter $
//
// 
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// G4MTcoutDestination.hh
//
// ---------------------------------------------------------------
#ifndef G4MTcoutDestination_H
#define G4MTcoutDestination_H

#include <iostream>
#include <sstream>
#include <fstream>

#include "globals.hh"
#include "G4MulticoutDestination.hh"
#include "G4StateManager.hh"

class G4LockcoutDestination;

class G4MTcoutDestination : public G4MulticoutDestination
{
  public:

    explicit G4MTcoutDestination(const G4int& threadId); 
    virtual ~G4MTcoutDestination(); 

    virtual void Reset(); 

    void SetDefaultOutput( G4bool addMasterDestination = true ,
                           G4bool formatAlsoMaster = true ); 

    void SetCoutFileName(const G4String& fileN = "G4cout.txt",
                               G4bool ifAppend = true); 
    void AddCoutFileName(const G4String& fileN = "G4cout.txt",
                               G4bool ifAppend = true);
    void SetCerrFileName(const G4String& fileN = "G4cerr.txt",
                               G4bool ifAppend = true); 
    void AddCerrFileName(const G4String& fileN = "G4cerr.txt",
                               G4bool ifAppend = true);

    void EnableBuffering(G4bool flag=true); 

    inline void SetPrefixString(const G4String& wd = "G4WT") { prefix = wd; } 

    void SetIgnoreCout(G4int tid = 0);
    inline void SetIgnoreInit(G4bool val=true) { ignoreInit = val; }

    inline G4String GetPrefixString() const { return prefix; } 
    inline G4String GetFullPrefixString() const
    { 
        std::stringstream os;
        os<<prefix<<id;
        return os.str();
    }

  protected:

    void AddMasterOutput( G4bool formatAlsoMaster); 
    void HandleFileCout( G4String fileN, G4bool appendFlag,
                         G4bool suppressDefault);
    void HandleFileCerr( G4String fileN, G4bool appendFlag,
                         G4bool suppressDefault);
  private:

    void DumpBuffer(); 
  
  private:

    // Reference to the default destination
    G4coutDestination* ref_defaultOut; 

    // Reference to the master destination
    G4coutDestination* ref_masterOut; 
    G4bool masterDestinationFlag;  
    G4bool masterDestinationFmtFlag; 

    const G4int id; 
    G4bool useBuffer; 
    G4bool ignoreCout; 
    G4bool ignoreInit; 

    G4String prefix;
    G4StateManager* stateMgr; 
};

#endif
