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
// G4WorkerRunManagerKernel
//
// Class description:
//
// This is a class for mandatory control of the Geant4 kernel.
// Implements Worker behavior in a MT application.
// This class is constructed by G4WorkerRunManager. If a user uses his/her
// own class instead of G4WorkerRunManager, this class must be instantiated
// at the very beginning of the application and must be deleted at the very
// end. Also, following methods must be invoked in the proper order:
//   DefineWorldVolume()
//   InitializePhysics()
//   RunInitialization()
//   RunTermination()
//
// User must provide his/her own classes derived from the following
// abstract class and register it to the RunManagerKernel:
//   G4VUserPhysicsList - Particle types, Processes and Cuts
//
// G4WorkerRunManagerKernel does not have any event loop. Handling of events
// is managed by G4RunManager.
//
// This class re-implements only the method that require special treatment
// to implement worker behavior

// Authors: M.Asai, A.Dotti (SLAC), 2013
// --------------------------------------------------------------------
#ifndef G4WorkerRunManagerKernel_hh
#define G4WorkerRunManagerKernel_hh 1

#include "G4RunManagerKernel.hh"

class G4WorkerRunManagerKernel : public G4RunManagerKernel
{
  public:

    G4WorkerRunManagerKernel();
    virtual ~G4WorkerRunManagerKernel();

  protected:

    void SetupShadowProcess() const;
      // Overwrite default behavior.
};

#endif
