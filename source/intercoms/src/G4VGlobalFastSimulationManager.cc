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
// $Id: G4VGlobalFastSimulationManager.cc,v 1.4 2002-11-20 14:46:01 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Abstract interface for GEANT4 Global Fast Simulation Manager.
// P. Mora de Freitas & M. Verderi 14/April/1999.

#include "G4VGlobalFastSimulationManager.hh"

G4VGlobalFastSimulationManager* 
G4VGlobalFastSimulationManager::fpConcreteInstance = 0;

G4VGlobalFastSimulationManager* 
G4VGlobalFastSimulationManager::GetConcreteInstance ()
{
  return fpConcreteInstance;
}

void
G4VGlobalFastSimulationManager::
SetConcreteInstance (G4VGlobalFastSimulationManager* m)
{
  fpConcreteInstance = m;
}
