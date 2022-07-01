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

#include "globals.hh"
#include "G4FastSimulationVector.hh"

#include "G4VGlobalFastSimulationManager.hh"
#include "G4FastSimulationManager.hh"
#include "G4FastSimulationManagerProcess.hh"


class G4FastSimulationMessenger;

enum  listType {
  NAMES_ONLY,
  MODELS,
  ISAPPLICABLE
};


// Class Description:
// This a singleton class which provides the management of the G4FastSimulationManager
// objects and some ghost facilities. 
//
// You can get access to it by:
//
// #include "G4GlobalFastSimulationManager.hh"
// ...
// ...
// G4GlobalFastSimulationManager* globalFSM;
// globalFSM = G4GlobalFastSimulationManager::getGlobalFastSimulationManager();
// ...
// ...
//    
// Presently, you will mainly need to use the GlobalFastSimulationManager if you use ghost 
// geometries.
//

class G4GlobalFastSimulationManager
{

public: // With  description 

  static G4GlobalFastSimulationManager* GetGlobalFastSimulationManager();
  // Provides a global access to the GlobalFastSimulationManager
  static G4GlobalFastSimulationManager* GetInstance();
  // Same as GetGlobalFastSimulationManager()
  
  G4VFastSimulationModel* GetFastSimulationModel(const G4String& modelName,
						 const G4VFastSimulationModel* previousFound = 0) const;
  // Iterative fetch of G4VFastSimulationModel objects by name:
  //    o returns the G4VFastSimulationModel* of model with name modelName;
  //    o returns 0 if no model found;
  //    o usage:
  //        myModel = gblManager->GetFastSimulationModel("MyModel");
  //    o note for the case of several models having the same name:
  //        - to get the first "MyModel" model:
  //             myModel1 = gblManager->GetFastSimulationModel("MyModel", 0);
  //        - to get the next one:
  //             myModel2 = gblManager->GetFastSimulationModel("MyModel", myModel1);
  //        - and so on.
  //        - When gblManager->GetFastSimulationModel("MyModel", myModel_n)
  //          returns a null pointer, no extra model with name "MyModel" exist.

   
public: // Without description

  // Destructor
  ~G4GlobalFastSimulationManager(); 

  //
  // G4FastSimulationManager(Process)'s management, no intended for general use.
  //
  // Methods for a G4FastSimulationManager to register itself
  //
  void    AddFastSimulationManager(G4FastSimulationManager*);
  void RemoveFastSimulationManager(G4FastSimulationManager*);
  //
  // G4FastSimulationManagerProcess bookeeping:
  //
  void    AddFSMP(G4FastSimulationManagerProcess*);
  void RemoveFSMP(G4FastSimulationManagerProcess*);


  // Flag that the Parameterisation must be closed.
  void FastSimulationNeedsToBeClosed();


public: // With  description 
  void ShowSetup();
  // Show the fast simulation setup : world(s), region(s), model(s) and links between them.
  // Requires the geometry to be closed.


public: // Without description

  void ListEnvelopes(const G4String&                 aName = "all",
		     listType                    aListType = NAMES_ONLY);
  void ListEnvelopes(const G4ParticleDefinition*                       );  
  
  void   ActivateFastSimulationModel(const G4String&);
  void InActivateFastSimulationModel(const G4String&);

  void Flush();

private:
  // Private construtor insures singleton class
  G4GlobalFastSimulationManager();

  // recursive display of regions, models, etc...
  void DisplayRegion(G4Region* motherRegion, G4int depth, std::vector<G4ParticleDefinition*>& particles) const;

  // The single instance.
  static G4ThreadLocal G4GlobalFastSimulationManager*           fGlobalFastSimulationManager;
  G4FastSimulationMessenger*                       fTheFastSimulationMessenger;
  G4FastSimulationVector <G4FastSimulationManager>             ManagedManagers;
  G4FastSimulationVector <G4FastSimulationManagerProcess>          fFSMPVector;
};

#endif 
// end of #ifndef G4GlobalFastSimulationManager_hh
