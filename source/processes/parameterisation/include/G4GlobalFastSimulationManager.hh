// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GlobalFastSimulationManager.hh,v 1.2 1999-04-15 17:17:20 mora Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
//  G4GlobalFastSimulationManager.hh
//
//  Description:
//    A singleton class which manages the Fast Simulation managers 
//    attached to envelopes.
//
//  History:
//    June 98: Verderi && MoraDeFreitas - "G4ParallelWorld" becomes
//             "G4FlavoredParallelWorld"; some method name changes;
//             GetFlavoredWorldForThis now returns a 
//             G4FlavoredParallelWorld pointer.
//    Feb 98: Verderi && MoraDeFreitas - First Implementation.
//
//---------------------------------------------------------------

#ifndef  G4GlobalFastSimulationManager_hh
#define  G4GlobalFastSimulationManager_hh

#include <rw/tpordvec.h>

#include "G4VGlobalFastSimulationManager.hh"
#include "G4FastSimulationManager.hh"
#include "G4FastSimulationManagerProcess.hh"
#include "G4StateManager.hh"
#include "G4VStateDependent.hh"
#include "G4FlavoredParallelWorld.hh"

enum  listType {
  NAMES_ONLY,
  MODELS,
  ISAPPLICABLE
};

class G4FastSimulationMessenger;

class G4GlobalFastSimulationManager : public G4VStateDependent, 
				      public G4VGlobalFastSimulationManager
{
public:
  // Global access to the GlobalFastSimulationManager
  static G4GlobalFastSimulationManager* GetGlobalFastSimulationManager();
  
  // Destructor
  ~G4GlobalFastSimulationManager(); 

  // G4FastSimulationManager's management :
  //
  // Methods for a G4FastSimulationManager to register itself
  //
  void AddFastSimulationManager(G4FastSimulationManager*);
  void RemoveFastSimulationManager(G4FastSimulationManager*);

  // Open/Close Parameterisation
  void FastSimulationNeedsToBeClosed();
  void CloseFastSimulation();

  // print/control commands
  void ListEnvelopes(const G4String& aName = "all",
		     listType aListType = NAMES_ONLY);
  void ListEnvelopes(const G4ParticleDefinition*);  
  
  void ActivateFastSimulationModel(const G4String&);
  void InActivateFastSimulationModel(const G4String&);

  // G4FastSimulationProcess interface
  G4FlavoredParallelWorld* GetFlavoredWorldForThis(G4ParticleDefinition *);

  // G4StateManager interface
  G4bool Notify(G4ApplicationState requestedState);

private:
  // Private construtor insures singleton class
  G4GlobalFastSimulationManager();

  // The single instance.
  static G4GlobalFastSimulationManager* fGlobalFastSimulationManager;

  // The G4FastSimulationMessenger
  G4FastSimulationMessenger* fTheFastSimulationMessenger;

  // List of G4FastSimulationManagers
  RWTPtrOrderedVector<G4FastSimulationManager> ManagedManagers;

  // fClosed flags if the NeededFlavoredWorlds List was Build.
  G4bool fClosed;

  // List of needed ParallelWorlds after close
  RWTPtrOrderedVector<G4FlavoredParallelWorld> NeededFlavoredWorlds;

  // Internal fonction to Build world volume clones.
  G4VPhysicalVolume* GiveMeAWorldVolumeClone();
};

#endif 
// end of #ifndef G4GlobalFastSimulationManager_hh
