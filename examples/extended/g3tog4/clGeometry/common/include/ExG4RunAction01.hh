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
//
/// \file ExG4RunAction01.hh
/// \brief Definition of the ExG4RunAction01 class

#ifndef ExG4RunAction01_h
#define ExG4RunAction01_h 1

#include "ExG4RunAction01Messenger.hh"

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

/// Run action class
///
/// It provides:
/// - setting verbose level
/// - option to save random number status at begin and end of run 
/// - option to read random number from a file at begin of run 
/// - option for automatic (time-based) random seeds for each run

class ExG4RunAction01 : public G4UserRunAction
{
  public:
    ExG4RunAction01();
    ~ExG4RunAction01();

  public:
    // methods
    virtual void BeginOfRunAction(const G4Run* run);
    virtual void EndOfRunAction(const G4Run* run);

    // set methods
    void  SetRndmFileName(const G4String& fileName);
    void  SetVerboseLevel(G4int level);
    void  SetSaveRndm(G4bool value);
    void  SetReadRndm(G4bool value);
    void  SetAutoSeed(G4bool value);

  private:
    // static data members
    static const G4String fgkDefaultRndmFileName;
  
    // data members
    ExG4RunAction01Messenger fMessenger;
    G4String fRndmFileName;    
    G4int    fVerboseLevel;
    G4bool   fSaveRndm;    
    G4bool   fReadRndm;
    G4bool   fAutoSeed;   
};
// inline functions

inline void ExG4RunAction01::SetVerboseLevel(G4int level) {
  fVerboseLevel = level;
}  

inline void ExG4RunAction01::SetRndmFileName(const G4String& fileName) {
  fRndmFileName = fileName;
}  

inline void ExG4RunAction01::SetSaveRndm(G4bool value) {
  fSaveRndm = value;
}

inline void ExG4RunAction01::SetReadRndm(G4bool value) {
  fReadRndm = value;
}

inline void ExG4RunAction01::SetAutoSeed(G4bool value) {
  fAutoSeed = value;
}

#endif
