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
// G4WrapperProcess
//
// Class description:
//
// Virtual class for wrapper process objects. 

// Author: H.Kurahige, 18 December 1996
// --------------------------------------------------------------------
#ifndef G4WrapperProcess_hh 
#define G4WrapperProcess_hh 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

class G4WrapperProcess : public G4VProcess
{
  public:

    G4WrapperProcess(const G4String& aName =  "Wrapped",
                     G4ProcessType aType = fNotDefined);
      // Constructor requires the process name and type

    G4WrapperProcess(const G4WrapperProcess& right);
      // Copy constructor copies the name but does not copy the 
      // physics table (null pointer is assigned)

    virtual ~G4WrapperProcess();
      // Destructor 

    G4WrapperProcess& operator=(const G4WrapperProcess&) = delete;

    inline G4bool operator==(const G4WrapperProcess& right) const;
    inline G4bool operator!=(const G4WrapperProcess& right) const;
      // Equality operators

    virtual void RegisterProcess(G4VProcess*);
    virtual const G4VProcess* GetRegisteredProcess() const;

    ////////////////////////////
    // DoIt    /////////////////
    ///////////////////////////

    virtual G4VParticleChange* PostStepDoIt(
                             const G4Track& track,
                             const G4Step&  stepData
                            );

    virtual G4VParticleChange* AlongStepDoIt(
                             const G4Track& track,
                             const G4Step& stepData
                            );
    virtual G4VParticleChange* AtRestDoIt(
                             const G4Track& track,
                             const G4Step& stepData
                            );

    //////////////////////////
    // GPIL    //////////////
    /////////////////////////  

    virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double  previousStepSize,
                             G4double  currentMinimumStep,
                             G4double& proposedSafety,
                             G4GPILSelection* selection);

    virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4ForceCondition* condition
                            );

    virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            );
  
    virtual G4bool IsApplicable(const G4ParticleDefinition&);
      // Returns true if this process object is applicable to
      // the particle type
      // Process will not be registered to a particle if IsApplicable is false   

    virtual void BuildPhysicsTable(const G4ParticleDefinition&);
      // Messaged by the Particle definition (via the Process manager)
      // whenever cross-section tables have to be rebuilt (i.e. if new
      // materials have been defined). 
      // It is overloaded by individual processes when they need physics
      // tables

    // Processes which build (for example in their constructors)
    // physics tables independent of cuts should preferably use a private
    // void BuildThePhysicsTable() function.
    // **Not** another BuildPhysicsTable().
 
    virtual void PreparePhysicsTable(const G4ParticleDefinition&);
      // Messaged by the Particle definition (via the Process manager)
      // whenever cross-section tables have to be prepare for rebuilt 
      // (i.e. if new materials have been defined). 
      // It is overloaded by individual processes when they need physics
      // tables


    virtual G4bool StorePhysicsTable(const G4ParticleDefinition* ,
                                     const G4String& directory, 
                                     G4bool ascii = false); 
      // Store PhysicsTable in a file. 
      // Return false in case of failure at I/O 
 
    virtual G4bool RetrievePhysicsTable( const G4ParticleDefinition* ,
                                         const G4String& directory, 
                                         G4bool ascii = false);
      // Retrieve Physics from a file. 
      // Return true if the Physics Table can be build by using file.
      // Return false if the process has no functionality or in case of failure.
      // File name should be defined by each process and the file should be
      // placed under the directory specified by the argument

    virtual void StartTracking(G4Track*);
    virtual void EndTracking();
      // inform Start/End of tracking for each track to the physics process 
 
    virtual void SetProcessManager(const G4ProcessManager*); 
      // A process manager sets its own pointer when the process is registered
      // in the process Manager
    virtual  const G4ProcessManager* GetProcessManager(); 
      // Get the process manager which the process belongs to
  
    virtual void ResetNumberOfInteractionLengthLeft();
      // Reset (determine the value of) NumberOfInteractionLengthLeft
    virtual void SetMasterProcess(G4VProcess* masterP);
      // Needed for MT, forward call to underlying process 

  protected:

    G4VProcess* pRegProcess = nullptr;
};

// ------------------------
// Inline operators
// ------------------------

inline
G4bool G4WrapperProcess::operator==(const G4WrapperProcess& right) const
{
  return (this == &right);
}

inline
G4bool G4WrapperProcess::operator!=(const G4WrapperProcess& right) const
{
  return (this != &right);
}

#endif
