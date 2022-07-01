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
// G4VGlobalFastSimulationManager
//
// Class description:
//
// Abstract interface for Global Fast Simulation Manager
// G4GlobalFastSimulationManager is a "Singleton".
// This class is an abstract interface for G4GlobalFastSimulationManager.
// It has the public access function GetConcreteInstance(), which is used
// to obtain a pointer to the concrete G4GlobalFastSimulationManager, should
// it exist. Then:
//
// G4VGlobalFastSimulationManager* pVFSMan =
//     G4VGlobalFastSimulationManager::GetConcreteInstance ();
//
// 'pVFSMan' points to the real (concrete) G4GlobalFastSimulationManager if
// at least a parameterisation envelope exists, otherwise is null.
//
// Thus all code must be protected, for example by:
//  if (pVFSMan)
//    G4FlavoredParallelWorld* = pVFSMan -> GetFlavoredWorldForThis(p);

// Authors: P. Mora de Freitas & M. Verderi, 14 April 1999
// --------------------------------------------------------------------
#ifndef G4VGLOBALFASTSIMULATIONMANAGER_HH
#define G4VGLOBALFASTSIMULATIONMANAGER_HH 1

#include "G4Types.hh"
#include "icomsdefs.hh"

class G4VFlavoredParallelWorld;
class G4ParticleDefinition;

class G4VGlobalFastSimulationManager
{
  public:

    static G4VGlobalFastSimulationManager* GetConcreteInstance();
      // Returns pointer to the actual Global Fast Simulation manager if
      // at least a parameterisation envelope exists

    virtual ~G4VGlobalFastSimulationManager() = default;

    virtual G4VFlavoredParallelWorld* GetFlavoredWorldForThis(
                                                     G4ParticleDefinition*) = 0;
      // VGlobalFastSimulationManager interface for visualisation

  protected:

    static void SetConcreteInstance(G4VGlobalFastSimulationManager*);
      // Sets the pointer to the actual Global Fast Simulation manager

    G4ICOMS_DLL
    static G4ThreadLocal G4VGlobalFastSimulationManager* fpConcreteInstance;
      // Pointer to real G4GlobalFastSimulationManager
};

#endif
