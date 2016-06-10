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
// $Id: $
//
// --------------------------------------------------------------------
// GEANT 4 class header file 
//
// Class Description:
//     A storage class for G4BiasingTrackData objects.
//
//      ----------------G4BiasingTrackDataStore ----------------
//
// Author: M.Verderi (LLR), November 2013
//
// ---------------------------------------------------------------------

#ifndef G4BiasingTrackDataStore_hh
#define G4BiasingTrackDataStore_hh

// A singleton that must be *thread local*

class G4BiasingTrackData;
class G4Track;
#include <map>
#include "G4ThreadLocalSingleton.hh"
class G4BiasingTrackDataStore {
  friend class G4ThreadLocalSingleton<G4BiasingTrackDataStore>;
public:
  static G4BiasingTrackDataStore* GetInstance();
  ~G4BiasingTrackDataStore();

  void   Register(G4BiasingTrackData*);
  void DeRegister(G4BiasingTrackData*);

  G4BiasingTrackData* GetBiasingTrackData( const G4Track* track ) {return fTrackDataStore[track]; }

  const std::map < const G4Track*, G4BiasingTrackData* >& GetMap() const {return fTrackDataStore;}

private:
  G4BiasingTrackDataStore();
  std::map < const G4Track*, G4BiasingTrackData* > fTrackDataStore;
};

#endif
