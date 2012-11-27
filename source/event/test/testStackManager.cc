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
// testStackManager.cc
//


#include "G4StackManager.hh"
#include "G4Trajectory.hh"
#include "G4Track.hh"
#include "G4ios.hh"

#include "time.h"

int main()
{
  //--- create a Stackmanager----------------------------------------
  G4StackManager stackmanager;
  G4VTrajectory* traj = new G4Trajectory();
  G4Track* track  = new G4Track(new G4DynamicParticle(), 0.0, G4ThreeVector());

  G4VTrajectory* rtraj;
  G4Track* rtrack;
  clock_t start, finish;

  start = clock();

  for(int i = 0; i < 10000; i++) stackmanager.PushOneTrack(track, traj);
  for(int i = 0; i < 100000000; i++) {
    rtrack = stackmanager.PopNextTrack(&rtraj);
    stackmanager.PushOneTrack(track, traj);
  }
  for(int i = 0; i < 10000; i++) rtrack = stackmanager.PopNextTrack(&rtraj);

  finish = clock();

  G4cout << "Execution time: " << double(finish - start)/CLOCKS_PER_SEC << " seconds" << G4endl;
  
  return 0;
}

