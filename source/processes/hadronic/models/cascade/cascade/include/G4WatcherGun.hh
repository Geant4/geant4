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
// $Id: G4WatcherGun.hh 66241 2012-12-13 18:34:42Z gunter $
//
// 20100407  M. Kelsey -- Return const-ref to avoid copy overhead.

#ifndef G4WATCHER_GUN_HH
#define G4WATCHER_GUN_HH

#include "G4NuclWatcher.hh"
#include <vector>

class G4WatcherGun {

public:

  G4WatcherGun();
  void setWatchers();

  const std::vector<G4NuclWatcher>& getWatchers() const { 
    return watchers; 
  };

private: 

  G4int verboseLevel;
  std::vector<G4NuclWatcher> watchers;

};        

#endif // G4WATCHER_GUN_HH 
