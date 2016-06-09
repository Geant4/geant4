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
// File: G4VHCIOentry.cc
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#include "G4VHCIOentry.hh"

// Addtional Include:
#include "G4HCIOcatalog.hh"

// Implementation of Constructor #1
G4VHCIOentry::G4VHCIOentry(std::string n)
 : m_name(n)
{
  G4HCIOcatalog* c = G4HCIOcatalog::GetHCIOcatalog();
  c->RegisterEntry(this);

  m_verbose = G4PersistencyCenter::GetPersistencyCenter()->VerboseLevel();
}

// End of G4VHCIOentry.cc

