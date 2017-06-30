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
// $Id: G4PDefManager.hh 103108 2017-03-16 13:00:35Z gcosmo $
//
//
// ------------------------------------------------------------
//
//	GEANT 4 class header file
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
#ifndef G4PDefManager_HH
#define G4PDefManager_HH

#include <stdlib.h>

#include "globals.hh"
#include "pwdefs.hh"
#include "G4AutoLock.hh"

class G4ProcessManager;

// G4PDefData is the private data from the object to be split
class G4PDefData
{
  // Encapsulates the fields of the class G4ParticleDefinition
  // that may not be read-only.

  public:

    void initialize();

    G4ProcessManager *theProcessManager;
};

// The type G4PDefManager is introduced to encapsulate the methods used by
// both the master thread and worker threads to allocate memory space for
// the fields encapsulated by the class G4PDefData. When each thread
// changes the value for these fields, it refers to them using a macro
// definition defined below. For every G4ParticleDefinition instance,
// there is a corresponding G4PDefData instance. All G4PDefData instances
// are organized by the class G4PDefManager as an array.
// The field "int g4particleDefinitionInstanceID" is added to the class G4ParticleDefinition.
// The value of this field in each G4ParticleDefinition instance is the
// subscript of the corresponding G4PDefData instance.
// In order to use the class G4PDefManager, we add a static member in the class
// G4ParticleDefinition as follows: "static G4PDefManager subInstanceManager".
// Both the master thread and worker threads change the length of the array
// for G4PDefData instances mutually along with G4ParticleDefinition
// instances are created. For each worker thread, it dynamically creates ions.
// Consider any thread A, if there is any other thread which creates an ion.
// This ion is shared by the thread A. So the thread A leaves an empty space
// in the array of G4PDefData instances for the ion.
//

class G4PDefManager
{
  public:

    G4PDefManager();
    G4int CreateSubInstance();
      // Invoked by the master or work thread to create a new subinstance
      // whenever a new split class instance is created. For each worker
      // thread, ions are created dynamically.


    void NewSubInstances();
      // Invoked by each worker thread to grow the subinstance array and
      // initialize each new subinstance using a particular method defined
      // by the subclass.

    void FreeSlave();
      // Invoked by all threads to free the subinstance array.

    G4PDefData*   GetOffset();

    void UseWorkArea( G4PDefData* newOffset ); // ,  G4int numObjects, G4int numSpace)

    G4PDefData* FreeWorkArea(); // G4int* numObjects, G4int* numSpace)

  public:

    G4PART_DLL G4ThreadLocalStatic G4int slavetotalspace;
    G4PART_DLL G4ThreadLocalStatic G4PDefData* offset;

  private:

    G4int totalobj;
    G4Mutex mutex;
};

#endif
