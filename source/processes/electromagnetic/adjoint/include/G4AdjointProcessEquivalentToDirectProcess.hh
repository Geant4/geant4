/////////////////////////////////////////////////////////////////////////////////
//      Class:		G4AdjointProcessEquivalentToDirectProcess
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	25 Sept. 2009 Created by L.Desorgher. Inspired from G4WrapperProcess  		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint process equivalent to direct process, used for some multiple scattering
//
// 


#ifndef G4AdjointProcessEquivalentToDirectProcess_h 
#define G4AdjointProcessEquivalentToDirectProcess_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

class G4AdjointProcessEquivalentToDirectProcess : public G4VProcess
{
  //  A virtual class for wrapper process objects.

  public: // with description
  //  constructor requires the process name and type
      G4AdjointProcessEquivalentToDirectProcess(const G4String& aName, G4VProcess* aProcess,G4ParticleDefinition* fwd_particle_def);

 

  public: 
  //  destructor 
      virtual ~G4AdjointProcessEquivalentToDirectProcess();

 
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
      virtual void      ResetNumberOfInteractionLengthLeft();
      // reset (determine the value of)NumberOfInteractionLengthLeft
   private:
      G4ParticleDefinition* theFwdParticleDef; 
      G4VProcess* theDirectProcess;  

};


#endif
