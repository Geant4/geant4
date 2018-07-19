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
#ifndef TBBTASK_HH
#define TBBTASK_HH

#include <tbb/task.h>
#include <tbb/concurrent_queue.h>
#include "G4Types.hh"

class G4Run;

class tbbTask : public tbb::task {
public:
  tbbTask(G4int anId, tbb::concurrent_queue<const G4Run*>* output=0 , G4int nEvts = 1 );
  virtual ~tbbTask();
  tbb::task* execute();
  
  unsigned int GetSlotId();  //
private:
  G4int m_nEvents;
  G4int m_taskID;
  tbb::concurrent_queue<const G4Run*>* m_output;
  G4bool m_beamOnCondition;
};

#endif //TBBTASK_HH