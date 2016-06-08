// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VGlobalFastSimulationManager.cc,v 1.1.8.1 1999/12/07 20:49:04 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// 
// Abstract interface for GEANT4 Global Fast Simulation Manager.
// P. Mora de Freitas & M. Verderi 14/April/1999.

#include "G4VGlobalFastSimulationManager.hh"

G4VGlobalFastSimulationManager* 
G4VGlobalFastSimulationManager::fpConcreteInstance = 0;
