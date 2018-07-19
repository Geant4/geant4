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
// $Id: G4FilecoutDestination.cc 103582 2017-04-18 17:24:45Z adotti $
//
// --------------------------------------------------------------------
//
// G4FilecoutDestination.cc
//
// Author: A.Dotti (SLAC), April 2017
// --------------------------------------------------------------------

#include "G4FilecoutDestination.hh"

#include <ios>

G4FilecoutDestination::~G4FilecoutDestination()
{
  Close();
  if ( m_output ) m_output.reset();
}

void G4FilecoutDestination::Open(std::ios_base::openmode mode)
{
  if ( m_name.isNull() )
  {
#ifndef __MIC
      // Cannot use G4Exception, because G4cout/G4cerr is not setup
      throw std::ios_base::failure("No output file name specified");
#endif
  }
  if ( m_output != nullptr && m_output->is_open() ) Close();
  m_output.reset( new std::ofstream(m_name , std::ios_base::out|mode) );
}

void G4FilecoutDestination::Close()
{
  if ( m_output && m_output->is_open() )
  {
    m_output->close();
  }
}

G4int G4FilecoutDestination::ReceiveG4cout(const G4String& msg)
{
  if ( m_output == nullptr || ! m_output->is_open() ) Open(m_mode);
  *m_output << msg;
  return 0;
}

G4int G4FilecoutDestination::ReceiveG4cerr(const G4String& msg)
{
  if ( m_output == nullptr || ! m_output->is_open() ) Open(m_mode);
  *m_output << msg;
  return 0;
}
