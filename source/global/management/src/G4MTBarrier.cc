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
// ---------------------------------------------------------------
/*
 * G4MTBarrier.cc
 *
 *  Created on: Feb 10, 2016
 *      Author: adotti
 *  Updated on: Feb 9, 2018
 *      Author: jmadsen
 */

#include "G4MTBarrier.hh"
#include "G4AutoLock.hh"

G4MTBarrier::G4MTBarrier(unsigned int numThreads ) :
  m_numActiveThreads(numThreads),
  m_counter(0)
{}

void G4MTBarrier::ThisWorkerReady() {
  //Step-1: Worker acquires lock on shared resource (the counter)
  G4AutoLock lock(&m_mutex);
  //Step-2: Worker increases counter
  ++m_counter;
  //Step-3: Worker broadcasts that the counter has changed
  G4CONDITIONBROADCAST(&m_counterChanged);
  //Step-4: Worker waits on condition to continue
  G4CONDITIONWAIT(&m_continue,&lock);
}

void G4MTBarrier::Wait() {
  while (true)
  {
      //Step-2: Acquires lock on shared resource (the counter)
      G4AutoLock lock(&m_mutex);
      //If the counter equals active threads, all threads are ready, exit the loop
      if ( m_counter == m_numActiveThreads ) { break; }
      //Step-3: Not all workers are ready, wait for the number to change
      //before repeating the check
      G4CONDITIONWAIT(&m_counterChanged,&lock);
  }
}

void G4MTBarrier::ReleaseBarrier() {
  //Step-4: re-aquire lock and re-set shared resource for future re-use
  G4AutoLock lock(&m_mutex);
  m_counter = 0;
  G4CONDITIONBROADCAST(&m_continue);
}

void G4MTBarrier::WaitForReadyWorkers() {
  //Step-1: Master enters a loop to wait all workers to be ready
  Wait();
  //Done, all workers are ready, broadcast a continue signal
  ReleaseBarrier();
}

void G4MTBarrier::ResetCounter() {
  G4AutoLock l(&m_mutex);
  m_counter = 0;
}

unsigned int G4MTBarrier::GetCounter() {
  G4AutoLock l(&m_mutex);
  const unsigned int result = m_counter;
  return result;
}
