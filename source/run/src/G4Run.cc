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
// $Id: G4Run.cc,v 1.9 2005/09/26 08:34:17 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

#include "G4Run.hh"

G4Run::G4Run()
:runID(0),numberOfEvent(0),numberOfEventToBeProcessed(0),HCtable(0),DCtable(0)
{;}

G4Run::~G4Run()
{;}

void G4Run::RecordEvent(const G4Event*)
{ numberOfEvent++; }

