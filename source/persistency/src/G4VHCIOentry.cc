//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
// File: G4VHCIOentry.cc
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#include "G4VHCIOentry.hh"

// Addtional Include:
#include "G4HCIOcatalog.hh"

// Implementation of Constructor #1
G4VHCIOentry::G4VHCIOentry(G4std::string n)
 : m_name(n)
{
  G4HCIOcatalog* c = G4HCIOcatalog::GetHCIOcatalog();
  c->RegisterEntry(this);

  m_verbose = G4PersistencyCenter::GetPersistencyCenter()->VerboseLevel();
}

// End of G4VHCIOentry.cc

