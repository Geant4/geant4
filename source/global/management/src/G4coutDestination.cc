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
//
//
// $Id: G4coutDestination.cc,v 1.1 2005/03/15 19:11:36 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// ----------------------------------------------------------------------
// G4coutDestination
//

#include "G4coutDestination.hh"

G4coutDestination::G4coutDestination()
{
}

G4coutDestination::~G4coutDestination()
{
}

G4int G4coutDestination::ReceiveG4cout(G4String)
{
  return 0;
}

G4int G4coutDestination::ReceiveG4cerr(G4String)
{
  return 0;
}
