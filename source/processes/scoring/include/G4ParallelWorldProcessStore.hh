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
// $Id: G4ParallelWorldProcessStore.hh 68733 2013-04-05 09:45:28Z gcosmo $
// GEANT4 tag $Name: geant4-09-04-ref-00 $
//
// 
//---------------------------------------------------------------
//
//  G4ParallelWorldProcessStore.hh
//
//  Description:
//    This procss takes a parallel world and limits a step
//   on the boundaries of volumes in the parallel world.
//    It invokes sensitive detectors assigned in the parallel
//   world.
//    It switches a material (and a region if defined) in the
//   assigned parallel world over the material (and the region)
//   in the mass world.
//
//---------------------------------------------------------------


#ifndef G4ParallelWorldProcessStore_h
#define G4ParallelWorldProcessStore_h 1

#include "globals.hh"
#include <map>
class G4ParallelWorldProcess;

//------------------------------------------
//
//        G4ParallelWorldProcessStore class
//
//------------------------------------------


// Class Description:

class G4ParallelWorldProcessStore 
: public std::map<G4ParallelWorldProcess*,G4String>
{
public: // with description
  static G4ParallelWorldProcessStore* GetInstance();
  static G4ParallelWorldProcessStore* GetInstanceIfExist();

private:
  static G4ThreadLocal G4ParallelWorldProcessStore* fInstance;

  //------------------------
  // Constructor/Destructor
  //------------------------
private:  
  G4ParallelWorldProcessStore();
public:
  virtual ~G4ParallelWorldProcessStore();
  
  //--------------------------------------------------------------
  // Set Paralle World
  //--------------------------------------------------------------

public:
  void SetParallelWorld(G4ParallelWorldProcess* proc,
                        G4String parallelWorldName);
  void UpdateWorlds();
  G4ParallelWorldProcess* GetProcess(G4String parallelWorldName);
  void Clear();

};

#endif
