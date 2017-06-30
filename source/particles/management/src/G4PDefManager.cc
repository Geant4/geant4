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
// $Id: G4PDefManager.cc 103108 2017-03-16 13:00:35Z gcosmo $
//
//
// ------------------------------------------------------------
//
//  GEANT4 class header file
//
// ---------------- G4PDefManager ----------------
//
// Utility template class for splitting RW data for thread-safety from
// classes: G4ParticleDefinition, G4VDecayChannel.
//
// ------------------------------------------------------------
// History:
// 01.25.2009 Xin Dong: First implementation from automatic MT conversion.
// ------------------------------------------------------------

#include <stdlib.h>

#include "G4PDefManager.hh"
#include "globals.hh"
#include "pwdefs.hh"
#include "G4AutoLock.hh"

void G4PDefData::initialize()
{
  theProcessManager = 0;
}

G4ThreadLocal G4int G4PDefManager::slavetotalspace=0;
G4ThreadLocal G4PDefData* G4PDefManager::offset=0;

G4PDefManager::G4PDefManager() : totalobj(0)
{
  G4MUTEXINIT(mutex);
}

G4int G4PDefManager::CreateSubInstance()
  // Invoked by the master or work thread to create a new subinstance
  // whenever a new split class instance is created. For each worker
  // thread, ions are created dynamically.
{
  G4AutoLock l(&mutex);
  totalobj++;
  if (totalobj > slavetotalspace)
  {
    l.unlock();
    NewSubInstances();
    l.lock();
  }
  return (totalobj - 1);
}

void G4PDefManager::NewSubInstances()
  // Invoked by each worker thread to grow the subinstance array and
  // initialize each new subinstance using a particular method defined
  // by the subclass.
{
  G4AutoLock l(&mutex);
  if (slavetotalspace  >= totalobj)  { return; }
  G4int originaltotalspace = slavetotalspace;
  slavetotalspace = totalobj + 512;
  offset = (G4PDefData *) realloc(offset, slavetotalspace * sizeof(G4PDefData));

  if (offset == 0)
  {
    G4Exception("G4PDefManager::NewSubInstances()",
                "OutOfMemory", FatalException, "Cannot malloc space!");
  }

  for (G4int i = originaltotalspace; i < slavetotalspace; i++)
  {
    offset[i].initialize();
  }
}

void G4PDefManager::FreeSlave()
  // Invoked by all threads to free the subinstance array.
{
  if (!offset)  { return; }
  free(offset);
  offset = 0;
}

G4PDefData* G4PDefManager::GetOffset()
{
  return offset;
}

void G4PDefManager::UseWorkArea( G4PDefData* newOffset )
  // Use recycled work area, which was created previously.
{
  if( offset && offset!=newOffset )
  {
    G4Exception("G4PDefManager::UseWorkspace()",
                "InvalidCondition", FatalException,
                "Thread already has workspace - cannot use another.");
  }
  offset= newOffset;
}

G4PDefData* G4PDefManager::FreeWorkArea()
  // Detach this thread from this Location
  // The object which calls this method is responsible for it.
{
  G4PDefData* offsetRet= offset;
  offset= 0;
  return offsetRet;
}
