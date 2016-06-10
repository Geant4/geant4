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
// $Id: G4WrapperProcess.hh 80787 2014-05-12 09:06:07Z gcosmo $
//
// 
// ------------------------------------------------------------
//        GEANT 4 class header file 
//
// Class Description
//
// This class is the virtual class for wrapper process objects. 

// ------------------------------------------------------------
//   New Physics scheme           18 Dec. 1996  H.Kurahige
// ------------------------------------------------------------

#ifndef G4WrapperProcess_h 
#define G4WrapperProcess_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

class G4WrapperProcess : public G4VProcess
{
  //  A virtual class for wrapper process objects.

  private:
  // hide default constructor and assignment operator as private 
  //  do not hide default constructor for alpha version 
      inline G4WrapperProcess & operator=(const G4WrapperProcess &right);

  public: // with description
  //  constructor requires the process name and type
      G4WrapperProcess(const G4String& aName =  "Wrapped",
                 G4ProcessType   aType = fNotDefined );

  //  copy constructor copys the name but does not copy the 
  //  physics table (0 pointer is assigned)
      G4WrapperProcess(const G4WrapperProcess &right);

  public: 
  //  destructor 
      virtual ~G4WrapperProcess();

  // equality opperators
      inline G4int operator==(const G4WrapperProcess &right) const;
      inline G4int operator!=(const G4WrapperProcess &right) const;

  public: // with description
    virtual void              RegisterProcess(G4VProcess*);
    virtual const G4VProcess* GetRegisteredProcess() const;

  protected:
    G4VProcess* pRegProcess;

  public: // with description
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
                            ) ;
  
  ////////////////////// 
      virtual G4bool IsApplicable(const G4ParticleDefinition&);
      // Returns true if this process object is applicable to
      // the particle type
      // Process will not be registered to a particle if IsApplicable is false   

      virtual void BuildPhysicsTable(const G4ParticleDefinition&);
      // Messaged by the Particle definition (via the Process manager)
      // whenever cross section tables have to be rebuilt (i.e. if new
      // materials have been defined). 
      // It is overloaded by individual processes when they need physics
      // tables. 

      // Processes which Build (for example in their
      // constructors) physics tables independent of cuts
      // should preferably use a
      // private void BuildThePhysicsTable()
      // function. Not another BuildPhysicsTable, please.
 
      virtual void PreparePhysicsTable(const G4ParticleDefinition&);
      // Messaged by the Particle definition (via the Process manager)
      // whenever cross section tables have to be prepare for rebuilt 
      // (i.e. if new materials have been defined). 
      // It is overloaded by individual processes when they need physics
      // tables. 

      // Processes which Build physics tables independent of cuts
      // (for example in their constructors)
      // should preferably use private 
      // void BuildThePhysicsTable() and void PreparePhysicsTable().
      // Not another BuildPhysicsTable, please.


      virtual G4bool StorePhysicsTable(const G4ParticleDefinition* ,
                                       const G4String& directory, 
                                       G4bool          ascii = false); 
      // Store PhysicsTable in a file. 
      // (return false in case of failure at I/O ) 
 
      virtual G4bool RetrievePhysicsTable( const G4ParticleDefinition* ,
                                           const G4String& directory, 
                                           G4bool          ascii = false);
      // Retrieve Physics from a file. 
      // (return true if the Physics Table can be build by using file)
      // (return false if the process has no functionality or in case of failure)
      // File name should be defined by each process 
      // and the file should be placed under the directory specifed by the argument. 
  ////////////////////////////
      virtual void StartTracking(G4Track*);
      virtual void EndTracking();
      // inform Start/End of tracking for each track to the physics process 
 
  public:
      virtual void SetProcessManager(const G4ProcessManager*); 
      // A process manager set its own pointer when the process is registered
      // the process Manager
      virtual  const G4ProcessManager* GetProcessManager(); 
      // Get the process manager which the process belongs to
  
   public:
      virtual void      ResetNumberOfInteractionLengthLeft();
      // reset (determine the value of)NumberOfInteractionLengthLeft
      virtual void SetMasterProcess(G4VProcess* masterP);
     // Needed for MT, forward call to underlying process 
};

inline
 G4WrapperProcess & G4WrapperProcess::operator=(const G4WrapperProcess &)
{
  G4Exception("G4WrapperProcess::operator=","Illegal operation",
              JustWarning,"Assignment operator is called");
  return *this;
}

inline
 G4int G4WrapperProcess::operator==(const G4WrapperProcess &right) const
{
  return (this == &right);
}

inline
 G4int G4WrapperProcess::operator!=(const G4WrapperProcess &right) const
{
  return (this !=  &right);
}

#endif
