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
/*
 * G4BuffercoutDestination.cc
 *
 *  Created on: Apr 14, 2017
 *      Author: adotti
 */

#include "G4BuffercoutDestination.hh"
#include "G4AutoLock.hh"
#include <iostream>


G4BuffercoutDestination::G4BuffercoutDestination(size_t max) :
  m_buffer_out("") , m_buffer_err(""), m_currentSize_out(0) ,
  m_currentSize_err(0), m_maxSize(max) {}

G4BuffercoutDestination::~G4BuffercoutDestination() {
  Finalize();
}

void G4BuffercoutDestination::Finalize() {
  FlushG4cerr();
  FlushG4cout();
}

G4int G4BuffercoutDestination::ReceiveG4cout(const G4String& msg) {
  m_currentSize_out += msg.size();
  m_buffer_out << msg;
  //If there is a max size and it has been reached, flush
  if ( m_maxSize>0 && m_currentSize_out >= m_maxSize ) {
      FlushG4cout();
  }
  return 0;
}

G4int G4BuffercoutDestination::ReceiveG4cerr(const G4String& msg) {
  m_currentSize_err += msg.size();
  m_buffer_err << msg;
  //If there is a max size and it has been reached, flush
  if ( m_maxSize>0 && m_currentSize_err >= m_maxSize ) {
      FlushG4cerr();
  }
  return 0;
}

G4int G4BuffercoutDestination::FlushG4cout() {
  std::cout<<m_buffer_out.str()<<std::flush;
  ResetCout();
  return 0;
}

void G4BuffercoutDestination::ResetCout() {
  m_buffer_out.str("");
  m_buffer_out.clear();
  m_currentSize_out = 0;
}

G4int G4BuffercoutDestination::FlushG4cerr() {
  std::cerr<<m_buffer_err.str()<<std::flush;
  ResetCerr();
  return 0;
}

void G4BuffercoutDestination::ResetCerr() {
  m_buffer_err.str("");
  m_buffer_err.clear();
  m_currentSize_err = 0;
}
