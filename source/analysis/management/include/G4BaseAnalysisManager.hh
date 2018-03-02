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
// $Id: G4BaseAnalysisManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

// Base class for the object managers.
// It handles first object ID and its lock and provides common utility methods.

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4BaseAnalysisManager_h
#define G4BaseAnalysisManager_h 1

#include "G4Fcn.hh"
#include "G4BinScheme.hh"
#include "globals.hh"

#include "G4AnalysisManagerState.hh"

class G4BaseAnalysisManager
{
  public:
    explicit G4BaseAnalysisManager(const G4AnalysisManagerState& state);
    virtual ~G4BaseAnalysisManager();
    
    // methods

    // The ids of objects are generated automatically
    // starting from 0; with the following function it is possible to
    // change the first Id to start from other value
    G4bool SetFirstId(G4int firstId);
    void  SetLockFirstId(G4bool lockFirstId);

    // Access method
    G4int GetFirstId() const;

  protected:
    // data members
    const G4AnalysisManagerState& fState;
    G4int    fFirstId;
    G4bool   fLockFirstId;
};

// inline functions

inline void G4BaseAnalysisManager::SetLockFirstId(G4bool lockFirstId) {
  fLockFirstId = lockFirstId;
}  

inline G4int G4BaseAnalysisManager::GetFirstId() const {
  return fFirstId;
}  

#endif

