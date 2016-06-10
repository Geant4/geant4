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
// 
// ------------------------------------------------------------
//
//	GEANT 4 class header file 
//
// ---------------- G4PDefSplitter ----------------
//
// Utility template class for splitting RW data for thread-safety from
// classes: G4ParticleDefinition, G4VDecayChannel.
//
// ------------------------------------------------------------
// History:
// 01.25.2009 Xin Dong: First implementation from automatic MT conversion.
// ------------------------------------------------------------
#ifndef G4PDEFSPLITTER_HH
#define G4PDEFSPLITTER_HH

#include <stdlib.h>

#include "globals.hh"
#include "pwdefs.hh"

template <class T>  // T is the private data from the object to be split
class G4PDefSplitter
{
  public:

    G4PDefSplitter() : totalobj(0) {}

    G4int CreateSubInstance()
      // Invoked by the master or work thread to create a new subinstance
      // whenever a new split class instance is created. For each worker
      // thread, ions are created dynamically.
    {
      totalobj++;
      if (totalobj > slavetotalspace)  { NewSubInstances(); }
      return (totalobj - 1);
    }

    void NewSubInstances()
      // Invoked by each worker thread to grow the subinstance array and
      // initialize each new subinstance using a particular method defined
      // by the subclass.
    {
      if (slavetotalspace  >= totalobj)  { return; }
      G4int originaltotalspace = slavetotalspace;
      slavetotalspace = totalobj + 512;
      offset = (T *) realloc(offset, slavetotalspace * sizeof(T));

      if (offset == 0)
      {
        G4Exception("G4PDefSplitter::NewSubInstances()",
                    "OutOfMemory", FatalException, "Cannot malloc space!");
      }

      for (G4int i = originaltotalspace; i < slavetotalspace; i++)
      {
        offset[i].initialize();
      }
    }

    void FreeSlave()
      // Invoked by all threads to free the subinstance array.
    {
      if (!offset)  { return; }
      free(offset);
      offset = 0;
    }

  public:

    G4PART_DLL static G4ThreadLocal G4int slavetotalspace;
    G4PART_DLL static G4ThreadLocal T* offset;

  private:

    G4int totalobj;
};

#endif
