//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4WrapperProcess.hh,v 1.1 2001-11-07 11:53:12 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
//
// Class Description
//  This class is the virtual class for wrapper process objects. 
//
// ------------------------------------------------------------
//   New Physics scheme           18 Dec. 1996  H.Kurahige
// ------------------------------------------------------------

#ifndef G4WrapperProcess_h 
#define G4WrapperProcess_h 1

#include "globals.hh"
#include "G4ios.hh"


#include "G4VProcess.hh"

class G4WrapperProcess  :public G4VProcess
{
  //  A virtual class for wrapper process objects.

  private:
  // hide default constructor and assignment operator as private 
  //  do not hide default constructor for alpha version 
      G4WrapperProcess & operator=(const G4WrapperProcess &right);

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

  // equal opperators
      G4int operator==(const G4WrapperProcess &right) const;
      G4int operator!=(const G4WrapperProcess &right) const;

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


      virtual G4bool StorePhysicsTable(G4ParticleDefinition* ,
				       const G4String& directory, 
				       G4bool          ascii = false); 
      // Store PhysicsTable in a file. 
      // (return false in case of failure at I/O ) 
 
      virtual G4bool RetrievePhysicsTable( G4ParticleDefinition* ,
					   const G4String& directory, 
				           G4bool          ascii = false);
      // Retrieve Physics from a file. 
      // (return true if the Physics Table can be build by using file)
      // (return false if the process has no functionality or in case of failure)
      // File name should be defined by each process 
      // and the file should be placed under the directory specifed by the argument. 
  ////////////////////////////
      virtual void StartTracking();
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

};


inline 
 void G4WrapperProcess::ResetNumberOfInteractionLengthLeft()
{
  pRegProcess->ResetNumberOfInteractionLengthLeft();
}


inline 
 G4double G4WrapperProcess::AlongStepGetPhysicalInteractionLength( const G4Track& track,
								   G4double  previousStepSize,
								   G4double  currentMinimumStep,
								   G4double& proposedSafety,
								   G4GPILSelection* selection     )
{
  return pRegProcess->AlongStepGetPhysicalInteractionLength( track,
							     previousStepSize,
							     currentMinimumStep,
							     proposedSafety,
							     selection     );
}

inline 
 G4double G4WrapperProcess::AtRestGetPhysicalInteractionLength( const G4Track& track,
								G4ForceCondition* condition )
{
  return pRegProcess->AtRestGetPhysicalInteractionLength( track,
				  condition );
}

inline 
 G4double G4WrapperProcess::PostStepGetPhysicalInteractionLength( const G4Track& track,
									G4double   previousStepSize,
									G4ForceCondition* condition )
{
   return pRegProcess->PostStepGetPhysicalInteractionLength( track,
							     previousStepSize,
							     condition );
}
      
inline 
 void G4WrapperProcess::SetProcessManager(const G4ProcessManager* procMan)
{
   pRegProcess->SetProcessManager(procMan); 
}

inline
 const G4ProcessManager* G4WrapperProcess::GetProcessManager()
{
  return     pRegProcess->GetProcessManager();
}

inline
 G4VParticleChange* G4WrapperProcess::PostStepDoIt(
						   const G4Track& track,
						   const G4Step&  stepData
						   )
{
  return     pRegProcess->PostStepDoIt( track, stepData );	
}

inline
 G4VParticleChange* G4WrapperProcess::AlongStepDoIt(
						    const G4Track& track,
						    const G4Step& stepData
						    )
{
 return     pRegProcess->AlongStepDoIt( track, stepData );	
}
 
inline
 G4VParticleChange* G4WrapperProcess::AtRestDoIt(
						 const G4Track& track,
						 const G4Step& stepData
						 )
{
 return     pRegProcess->AtRestDoIt( track, stepData );	
}

inline
 G4bool G4WrapperProcess::IsApplicable(const G4ParticleDefinition& particle)
{
  return     pRegProcess->IsApplicable(particle);
}

inline
 void G4WrapperProcess::BuildPhysicsTable(const G4ParticleDefinition& particle)
{
  return     pRegProcess->BuildPhysicsTable(particle);
}

inline
 G4bool G4WrapperProcess::StorePhysicsTable(G4ParticleDefinition* particle,
				       const G4String& directory, 
				       G4bool          ascii)
{
  return pRegProcess->StorePhysicsTable(particle,  directory,  ascii);
} 
 
inline
 G4bool G4WrapperProcess::RetrievePhysicsTable( G4ParticleDefinition* particle,
				       const G4String& directory, 
				       G4bool          ascii)
{
  return pRegProcess->RetrievePhysicsTable(particle,  directory,  ascii);
}  

inline
 void G4WrapperProcess::StartTracking()
{
  pRegProcess->StartTracking();
}

inline
 void G4WrapperProcess::EndTracking()
{
  pRegProcess->EndTracking();
}

inline
 void   G4WrapperProcess::RegisterProcess(G4VProcess* process)
{
  pRegProcess=process;
  theProcessName += process->GetProcessName();
  theProcessType = process->GetProcessType();
}

inline
 const G4VProcess* G4WrapperProcess::GetRegisteredProcess() const
{
  return pRegProcess;
} 

inline
 G4WrapperProcess::G4WrapperProcess(const G4String& aName,
				    G4ProcessType   aType)
  : G4VProcess(aName,aType), pRegProcess((G4VProcess*)(0))
{
}

inline
G4WrapperProcess::G4WrapperProcess(const G4WrapperProcess& right)
  : G4VProcess(*((G4VProcess*)(&right))), pRegProcess(right.pRegProcess)
{
}

inline
 G4WrapperProcess::~G4WrapperProcess()
{
  if (pRegProcess!=0) delete pRegProcess;
}

inline
 G4WrapperProcess & G4WrapperProcess::operator=(const G4WrapperProcess &)
{
  G4Exception("G4WrapperProcess::assignment operator is called");
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




















