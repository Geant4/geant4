// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VGlobalFastSimulationManager.cc,v 1.1 1999-04-15 15:43:19 mora Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Abstract interface for GEANT4 Global Fast Simulation Manager.
// P. Mora de Freitas & M. Verderi 14/April/1999.

#include "G4VGlobalFastSimulationManager.hh"

G4VGlobalFastSimulationManager* 
G4VGlobalFastSimulationManager::fpConcreteInstance = 0;
