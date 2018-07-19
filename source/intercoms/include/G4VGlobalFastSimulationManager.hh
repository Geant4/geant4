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
//
// $Id: G4VGlobalFastSimulationManager.hh 106172 2017-09-15 13:03:57Z gcosmo $
//
// 
// Abstract interface for GEANT4 Global Fast Simulation Manager.
// P. Mora de Freitas & M. Verderi 14/April/1999.
//
// Class description:
//
// G4GlobalFastSimulationManager is a "Singleton", i.e., only one instance 
// of it may exist. This is ensured by making the constructor private.
//
// G4VGlobalFastSimulationManager is an abstract interface for the
// G4GlobalFastSimulationManager one. It has the public access function
// GetConcreteInstance(), which is used to obtain a pointer to the concrete 
// G4GlobalFastSimulationManager, should it exist. After
//
// G4VGlobalFastSimulationManager* pVFSMan =  
//     G4VGlobalFastSimulationManager::GetConcreteInstance ();
//
// pVFSMan points to the real (concrete) G4GlobalFastSimulationManager if
// at least a parameterisation envelope exists, otherwise is zero.  
//
// Thus all code must be protected, for example by:
//   if (pVFSMan) 
//    G4FlavoredParallelWorld* =
//      pVFSMan -> GetFlavoredWorldForThis(p);
//

#ifndef G4VGLOBALFASTSIMULATIONMANAGER_HH
#define G4VGLOBALFASTSIMULATIONMANAGER_HH

#include "G4Types.hh"
#include "icomsdefs.hh"

class G4VFlavoredParallelWorld;
class G4ParticleDefinition;

class G4VGlobalFastSimulationManager
{

public:  // with description

  static G4VGlobalFastSimulationManager* GetConcreteInstance ();
    // Returns pointer to actual Global Fast Simulation manager if
    // at least a parameterisation envelope exists. Always check value.

  virtual ~G4VGlobalFastSimulationManager () {}

  virtual
  G4VFlavoredParallelWorld* GetFlavoredWorldForThis(G4ParticleDefinition*)=0;
    // VGlobalFastSimulationManager interface for visualisation.

protected:

  static void SetConcreteInstance (G4VGlobalFastSimulationManager*);
    // Sets the pointer to actual Global Fast Simulation manager.

  G4ICOMS_DLL static G4ThreadLocal G4VGlobalFastSimulationManager* fpConcreteInstance;  
    // Pointer to real G4GlobalFastSimulationManager.

};

#endif
