// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VUserPhysicsList.hh,v 1.7 2000-11-08 10:01:59 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
// Class Description:
//      This class is an abstruct class for
//      constructing particles and processes.
//      User must implement following four virtual methods
//      in his own concrete class derived from this class. 
//        G4VUserPhysicsList::ConstructParticle() 
//           Construct particles
//        G4VUserPhysicsList::constructPhysics() 
//           Construct procesess and register them to particles
//        G4VUserPhysicsList::SetCuts()
//           set cut values in range to all particles
//           (and rebuilding physics table will be invoked )
//
// ------------------------------------------- 
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
//          removed ConstructAllParticles() and related methods  
//          changed SetCuts method argument
//       modified                           08, Nov 2000 by H.Kurashige
//          added   Retrieve/StorePhysicsTable and related methods
// ------------------------------------------------------------
#ifndef G4VUserPhysicsList_h
#define G4VUserPhysicsList_h 1
#include "globals.hh"
#include "G4ios.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh" 

class G4UserPhysicsListMessenger;
class G4VProcess;

class G4VUserPhysicsList
{
  public: 
    G4VUserPhysicsList();
    virtual ~G4VUserPhysicsList();

  public:  // with description
    // By calling the "Construct" method, 
    // particles and processes are created    
    void Construct();
 

  protected: // with description
    // These two methods of  ConstructParticle() and ConstructProcess()
    // will be invoked in the Construct() method. 

    // each particle type will be instantiated
    virtual void ConstructParticle() = 0;
 
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess() = 0;

  protected: // with description
   //  User must invoke this method in his ConstructProcess() 
   //  implementation in order to insures particle transportation.
   //  !! Caution: this class must not be overriden !!
    void AddTransportation();

  /////////////////////////////////////
  public:  
   //  "SetCuts" method sets a cut value for all particle types 
   //   in the particle table
   virtual void SetCuts() = 0; 

  protected: // with description
   //  "SetCutsWithDefault" method sets a cut value with the default
   //   cut values for all particle types in the particle table
   void SetCutsWithDefault();   

  public: // with description
    // Invoke BuildPhysicsTable for all processes for each particle
    // In case of "Retrieve" flag is ON, PhysicsTable will be
    // retrieved from files
    void BuildPhysicsTable(G4ParticleDefinition* );    
 
    // Store PhysicsTable together with both material and cut value 
    // information in files under the specified directory.
    //  (return true if files are sucessfully created)
    G4bool  StorePhysicsTable(const G4String& directory = ".");

    // Return true if "Retrieve" flag is ON. 
    // (i.e. PhysicsTable will be retrieved from files)
    G4bool  IsPhysicsTableRetrieved() const;
    // Get directory path for physics table files.
    const G4String& GetPhysicsTableDirectory() const;

    // Set "Retrieve" flag
    // Directory path can be set together.
    // Null string (default) means directory is not changed 
    // from the current value 
    void    SetPhysicsTableRetrieved(const G4String& directory = "");
   
    // Reset "Retrieve" flag
    void    ResetPhysicsTableRetrieved();
  
  protected:  // with description
    // do BuildPhysicsTable for make the integral schema
    void BuildIntegralPhysicsTable(G4VProcess* ,G4ParticleDefinition*  );   

  protected: // with description
    // Retrieve PhysicsTable from files for proccess belongng the particle.
    // Normal BuildPhysics procedure of processes will be invoked, 
    // if it fails (in case of Process's RetrievePhysicsTable returns false)
    virtual void  RetrievePhysicsTable(G4ParticleDefinition* );

    // Store material information in files under the specified directory.
    virtual G4bool  StoreMaterialInfo(const G4String& directory);
    // Store cut values information in files under the specified directory.
    virtual G4bool  StoreCutValues(const G4String& directory);

    // check stored material and cut values
    virtual G4bool CheckForRetrievePhysicsTable(const G4String& directory);
    // check stored material is consistent with the current detector setup. 
    virtual G4bool  CheckMaterialInfo(const G4String& directory);
    // check stored cuvalue is consistent with the current detector setup. 
    virtual G4bool  CheckCutValues(const G4String& directory);

  protected: // with description
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

  public:  // with description
    //  set/get the default cut value
    //  Calling SetDefaultCutValue causes re-calcuration of cut values
    //  and physics tables just before the next event loop
    void     SetDefaultCutValue(G4double newCutValue);
    G4double GetDefaultCutValue() const;

  protected:
    // this is the default cut value for all particles
    G4double defaultCutValue;

  /////////////////////////////////////
  public: // with description
    // Print out the List of registered particles types
    void DumpList() const;

  public: // with description
  // Print out information of cut values
    void DumpCutValuesTable() const;
    void DumpCutValues(const G4String &particle_name = "ALL") const;
    void DumpCutValues(G4ParticleDefinition* ) const;

  protected: // with description
    // adds new ProcessManager to all particles in the Particle Table
    //   this routine is used in Construct()
    void InitializeProcessManager();

  public: // with description
    // remove and delete ProcessManagers for all particles in tha Particle Table
    //    this routine is invoked from RunManager 
    void RemoveProcessManager();

  public: // with description
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

 public: // with description
   void  SetVerboseLevel(G4int value);
   G4int GetVerboseLevel() const;
   // set/get controle flag for output message
   //  0: Silent
   //  1: Warning message
   //  2: More

 protected:
   G4int verboseLevel;

 protected:
   // flag to determine physics table will be build from file or not
   G4bool fRetrievePhysicsTable;  
 
   G4bool fIsCheckedForRetrievePhysicsTable;
   
   // directory name for physics table files 
   G4String directoryPhysicsTable;   

   // number of materials in G4MaterialTable
   // (this member is used by store/restore physics table)
   G4int numberOfMaterial;   
};


inline void G4VUserPhysicsList::Construct()
{
  if (verboseLevel >1) G4cout << "G4VUserPhysicsList::Construct()" << G4endl;  

  if (verboseLevel >1) G4cout << "Construct particles " << G4endl;  
  ConstructParticle();

  InitializeProcessManager();

  if (verboseLevel >1) G4cout << "Construct processes " << G4endl;  
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
    G4cout << " Verbose level is set to " << verboseLevel << G4endl;
  }
}

inline  G4int G4VUserPhysicsList::GetVerboseLevel() const
{
  return  verboseLevel;
}

inline  
 G4bool  G4VUserPhysicsList::IsPhysicsTableRetrieved() const
{
  return fRetrievePhysicsTable;  
}

inline 
  const G4String& G4VUserPhysicsList::GetPhysicsTableDirectory() const
{
  return directoryPhysicsTable;  
}

inline 
void  G4VUserPhysicsList::SetPhysicsTableRetrieved(const G4String& directory)
{
  fRetrievePhysicsTable = true;
  if(!directory.isNull()) {
    directoryPhysicsTable = directory;
  }
  fIsCheckedForRetrievePhysicsTable=false;
}
   
inline 
 void  G4VUserPhysicsList::ResetPhysicsTableRetrieved()
{
  fRetrievePhysicsTable = false;
  fIsCheckedForRetrievePhysicsTable=false;
}
#endif







