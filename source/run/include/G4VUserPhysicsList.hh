// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VUserPhysicsList.hh,v 1.3 1999-04-14 10:45:56 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This class is an abstruct class for
//      constructing particles and processes.
//      User must implement following four virtual methods
//      in his own concrete class derived from this class. 
//        G4VUserPhysicsList::ConstructParticle() 
//           Construct particles
//        G4VUserPhysicsList::constructPhysics() 
//           Construct procesess and register them to particles
//        G4VUserPhysicsList::SetCuts(G4double aValue)
//           set a cut value in range to all particles
//           (and rebuilding physics table will be invoked)
// 
//	History
//        first version                   09 Jan. 1998 by H.Kurashige 
//        modified                        24 Jan. 1998 by H.Kurashige
//          rename DumpCutValues/DumpCutValuesTable 
//          change SetCuts method
//          add    SetCutsWithDefault method
//        modified                       06 June 1998 by H.Kurashige
//          add    AddProcessManager
//          add    BuildPhysicsTable
//       modified                        29  June 1998 by H.Kurashige
//          add    AddProcessManager
//       modified                        05 Dec. 1998 by H.Kurashige
//          add    ConstructAllParticles()
//        modified                        14, Apr 1999 by H.Kurashige
//          change BuildPhysicsTable as public
// ------------------------------------------------------------
#ifndef G4VUserPhysicsList_h
#define G4VUserPhysicsList_h 1
#include "globals.hh"
#include "G4ios.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh" 

class G4UserPhysicsListMessenger;

class G4VUserPhysicsList
{
  public:
    G4VUserPhysicsList();
    virtual ~G4VUserPhysicsList();

  public:
    // By calling the "Construct" method, 
    // particles and processes are created    
    void Construct();
 
   protected:
    // These two methods of  ConstructParticle() and ConstructProcess()
    // will be invoked in the Construct() method. 

    // each particle type will be instantiated
    virtual void ConstructParticle() = 0;
 
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess() = 0;

  protected:
   //  User must invoke this method in his ConstructProcess() 
   //  implementation in order to insures particle transportation.
   //  !! Caution: this class must not be overriden !!
    void AddTransportation();

   // Construct All particles 
    virtual void ConstructAllParticles();
    // these methods Construct all particles in each category
    virtual void ConstructAllBosons();
    virtual void ConstructAllLeptons();
    virtual void ConstructAllMesons();
    virtual void ConstructAllBarions();
    virtual void ConstructAllIons();
    virtual void ConstructAllShortLiveds();

  public:  
   //  "SetCuts" method sets a cut value for all particle types 
   //   in the particle table
    virtual void SetCuts(G4double aCut) = 0;

   //  "SetCutsWithDefault" method sets a cut value with the default
   //   cut values for all particle types in the particle table
    void SetCutsWithDefault();   

  public:
    // 
    void BuildPhysicsTable(G4ParticleDefinition* );    
 
  protected:
    // Following are utility methods for SetCuts/reCalcCuts  

    // Reset cut values in energy for all particle types
    // By calling this methods, the run manager will invoke
    // SetCuts() just before event loop  
    void ResetCuts();

    // SetCutValue sets a cut value for a particle type
    void SetCutValue(G4double aCut, const G4String& name); 
    void ReCalcCutValue(const G4String& name); 

    //  "setCutsForOthers" method sets a cut value to all particle types 
    //  which have not be called SetCuts() methods yet.
    //  (i.e. particles which have no definit cut values)
    void SetCutValueForOthers(G4double cutValue);

    // "setCutsForOtherThan"  sets a cut value to all particle types
    // other than particle types specified in arguments
    void SetCutValueForOtherThan(G4double cutValue,
                                 G4ParticleDefinition* first,
                                 G4ParticleDefinition* second  = NULL,
				 G4ParticleDefinition* third   = NULL,
                                 G4ParticleDefinition* fourth  = NULL,
				 G4ParticleDefinition* fifth   = NULL,
                                 G4ParticleDefinition* sixth   = NULL,
				 G4ParticleDefinition* seventh = NULL,
                                 G4ParticleDefinition* eighth  = NULL,
				 G4ParticleDefinition* nineth  = NULL,
                                 G4ParticleDefinition* tenth   = NULL  );

    //  "reCalcCutsForOthers" method re-calculates a cut value 
    //  to all particle types which have not be called SetCuts() methods yet.
    void ReCalcCutValueForOthers();

  public: 
    //  set/get the default cut value
    //  Calling SetDefaultCutValue causes re-calcuration of cut values
    //  and physics tables just before the next event loop
    void     SetDefaultCutValue(G4double newCutValue);
    G4double GetDefaultCutValue() const;

  protected:
    // this is the default cut value for all particles
    G4double defaultCutValue;

  public:
    // Print out the List of registered particles types
    void DumpList() const;

  public:
  // Print out information of cut values
    void DumpCutValuesTable() const;
    void DumpCutValues(const G4String &particle_name = "ALL") const;
    void DumpCutValues(G4ParticleDefinition* ) const;

  protected:
    // adds new ProcessManager to all particles in the Particle Table
    //   this routine is used in Construct()
    void InitializeProcessManager();

  public:
    // remove and delete ProcessManagers for all particles in tha Particle Table
    //    this routine is invoked from RunManager 
    void RemoveProcessManager();

  public:
    // add process manager for particles created on-the-fly 
    void AddProcessManager(G4ParticleDefinition* newParticle,
			   G4ProcessManager*    newManager = NULL );
 
  protected:
    // the particle table has the complete List of existing particle types
    G4ParticleTable* theParticleTable;
    G4ParticleTable::G4PTblDicIterator* theParticleIterator;

  protected:
  // pointer to G4UserPhysicsListMessenger
    G4UserPhysicsListMessenger* theMessenger;

 public:
   void  SetVerboseLevel(G4int value);
   G4int GetVerboseLevel() const;

 protected:
   G4int verboseLevel;
   // controle flag for output message
   //  0: Silent
   //  1: Warning message
   //  2: More

};

inline void G4VUserPhysicsList::SetCutsWithDefault()
{
  SetCuts(defaultCutValue);
}  

inline void G4VUserPhysicsList::Construct()
{
  if (verboseLevel >1) G4cout << "G4VUserPhysicsList::Construct()" << endl;  

  if (verboseLevel >1) G4cout << "Construct particles " << endl;  
  ConstructParticle();

  InitializeProcessManager();

  if (verboseLevel >1) G4cout << "Construct processes " << endl;  
  ConstructProcess();
}

inline G4double G4VUserPhysicsList::GetDefaultCutValue() const
{
  return defaultCutValue;
}

inline void G4VUserPhysicsList::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
  if (verboseLevel >1){
    G4cout << "G4VUserPhysicsList::SetVerboseLevel  :";
    G4cout << " Verbose level is set to " << verboseLevel << endl;
  }
}

inline  G4int G4VUserPhysicsList::GetVerboseLevel() const
{
  return  verboseLevel;
}

#endif
