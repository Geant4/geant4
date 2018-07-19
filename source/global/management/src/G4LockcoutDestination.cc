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
// $Id: G4LockcoutDestination.cc 103582 2017-04-18 17:24:45Z adotti $
//
// --------------------------------------------------------------------
//
// G4LockcoutDestination.cc
//
// Author: A.Dotti (SLAC), April 2017
// --------------------------------------------------------------------

#include "G4LockcoutDestination.hh"
#include "G4AutoLock.hh"

namespace
{
  G4Mutex out_mutex = G4MUTEX_INITIALIZER;
}

G4LockcoutDestination::~G4LockcoutDestination()
{
}

G4int G4LockcoutDestination::ReceiveG4cout(const G4String& msg )
{
  G4AutoLock l(&out_mutex);
  // Forward call to base class
  return G4coutDestination::ReceiveG4cout(msg);
}

G4int G4LockcoutDestination::ReceiveG4cerr(const G4String& msg )
{
  G4AutoLock l(&out_mutex);
  return G4coutDestination::ReceiveG4cerr(msg);
}
