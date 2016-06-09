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
/// \file ExG4EventAction01.hh
/// \brief Definition of the ExG4EventAction01 class

#ifndef ExG4EventAction01_h
#define ExG4EventAction01_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "ExG4EventAction01Messenger.hh"

/// Event action class
///
/// It provides:
/// - setting verbose level and event modulo pronting
/// - option to save random number status at the end of event 

class ExG4EventAction01 : public G4UserEventAction
{
  public:
    ExG4EventAction01();
    virtual ~ExG4EventAction01();

  public:
    // methods from base class
    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    // set methods
    void SetVerboseLevel(G4int level);
    void SetPrintModulo(G4int value);
    void SetSaveRndm(G4bool value);
    void SetDrawFlag(G4String value);
        
  private:
    // static data members (constants)
    static const G4int  fgkDefaultVerboseLevel; 
    static const G4int  fgkDefaultPrintModulo; 
    
    // data members
    ExG4EventAction01Messenger  fMessenger;
    G4int     fVerboseLevel;
    G4int     fPrintModulo; 
    G4bool    fSaveRndm;  
};

// inline functions

inline void ExG4EventAction01::SetVerboseLevel(G4int level) {
  fVerboseLevel = level;
}  

inline void ExG4EventAction01::SetPrintModulo(G4int value) {
  fPrintModulo = value;
}

inline void ExG4EventAction01::SetSaveRndm(G4bool value) {
  fSaveRndm = value;
}

#endif
