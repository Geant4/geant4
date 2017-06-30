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
// $Id: G4VUserPhysicsList.hh 103803 2017-04-27 14:03:05Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
// Class Description:
//      This class is an abstract class for
//      constructing particles and processes.
//      User must implement following three virtual methods
//      in his/her own concrete class derived from this class. 
//        G4VUserPhysicsList::ConstructParticle() 
//           Construct particles
//        G4VUserPhysicsList::ConstructProcess() 
//           Construct procesess and register them to particles
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
//       modified                           08, Mar 2001 by H.Kurashige
//          added   binary mode for Retrieve/StorePhysicsTable
//          added   RetrieveCutValues and related
//          added   Set/ResetStoredInAscii() to switch on ascii mode 
//                  for Retrieve/StorePhysicsTable
//       modified for CUTS per REGION      10, Oct 2002 by H.Kurashige 
//          removed following methods 
//           void ReCalcCutValue() 
//           void SetCutValueForOthers()
//           void SetCutValueForOtherThan()	
//           void ReCalcCutValueForOthers()
//           virtual G4bool  StoreMaterialInfo()
//           virtual G4bool  StoreCutValues()
//           virtual G4bool  RetrieveCutValues()
//           virtual G4bool  CheckForRetrievePhysicsTable()
//           virtual G4bool  CheckMaterialInfo()
//          added    void BuildPhysicsTable()    
//       Added PhysicsListHelper           29 Apr. 2011 H.Kurashige
//       Added default impelmentation of SetCuts 10 June 2011 H.Kurashige 
//           SetCuts is not 'pure virtual' any more
//       Trasnformations for multi-threading 26 Mar. 2013 A. Dotti
//	 Added destructions 21 Apr 2017 A. Dotti
// ------------------------------------------------------------
#ifndef G4VUserPhysicsList_h
#define G4VUserPhysicsList_h 1

#include "globals.hh"
#include "tls.hh"
#include "rundefs.hh"
#include "G4ios.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh" 
#include "G4ProductionCutsTable.hh"
#include "G4VUPLSplitter.hh"

#include "G4Threading.hh"

class G4UserPhysicsListMessenger;
class G4PhysicsListHelper;
class G4VProcess;

class G4VUPLData
{
    //Encapsulate the fields of class G4VUserPhysicsList
    //that are per-thread.
public:
    void initialize();
    G4ParticleTable::G4PTblDicIterator* _theParticleIterator;
    G4UserPhysicsListMessenger* _theMessenger;
    G4PhysicsListHelper* _thePLHelper;
    G4bool _fIsPhysicsTableBuilt;
    G4int _fDisplayThreshold;
};
// The type G4VUPLManager is introduced to encapsulate the methods used by
// both the master thread and worker threads to allocate memory space for
// the fields encapsulated by the class G4VUPLData. When each thread
// changes the value for these fields, it refers to them using a macro
// definition defined below. For every G4VUserPhysicsList instance,
// there is a corresponding G4VUPLData instance. All G4VUPLData instances
// are organized by the class G4VUPLManager as an array.
// The field "int g4vuplInstanceID" is added to the class G4VUserPhysicsList.
// The value of this field in each G4VUserPhysicsList instance is the
// subscript of the corresponding G44VUPLData instance.
// In order to use the class G44VUPLManager, we add a static member in the class
// G4VUserPhysicsList as follows: "static G4VUPLManager subInstanceManager".
// Both the master thread and worker threads change the length of the array
// for G44VUPLData instances mutually along with G4VUserPhysicsList
// instances are created. For each worker thread, it dynamically creates ions.
// Consider any thread A, if there is any other thread which creates an ion.
// This ion is shared by the thread A. So the thread A leaves an empty space
// in the array of G4PDefData instances for the ion.
//
// Important Note: you may wonder why we are introducing this mechanism
//                 since there is only one PL for each application.
//                 This is true, in the sense that only one PL is allowed
//                 to be associated to a G4RunManager, however user can instantiate
//                 as many PLs are needed and at run-time select one of the PLs to be used
//                 we thus need this mechanism to guarantee that the system works without
//                 problems in case of this (unusual) case. This may be reviewed in the future
typedef G4VUPLSplitter<G4VUPLData> G4VUPLManager;
typedef G4VUPLManager G4VUserPhysicsListSubInstanceManager;

// This macros change the references to fields that are now encapsulated
// in the class G4VUPLData.
//
// Note1:  the use of this-> this is needed to avoid compilation errors
// when using templated class with T=G4VUserPhysicsList. Don't know why.
// Note2: the name of the first #define is different, because otherwise
//        we need to change its use in all classes that inherits from
//        this base class (all examples). However one should note comment
//        on JIRA task: http://jira-geant4.kek.jp/browse/DEV-27
//#define theParticleIterator ((this->subInstanceManager.offset[this->g4vuplInstanceID])._theParticleIterator)

class G4VUserPhysicsList
{
  public: 
    G4VUserPhysicsList();
    virtual ~G4VUserPhysicsList();

  // copy constructor and assignment operator
    G4VUserPhysicsList(const G4VUserPhysicsList&);
    G4VUserPhysicsList & operator=(const G4VUserPhysicsList&);

  public:  // with description
   // Each particle type will be instantiated
   // This method is invoked by the RunManger 
   virtual void ConstructParticle() = 0;

   // By calling the "Construct" method, 
   // process manager and processes are created. 
   void Construct();
 
   // Each physics process will be instantiated and
   // registered to the process manager of each particle type 
   // This method is invoked in Construct method 
   virtual void ConstructProcess() = 0;

  protected: // with description
   //  User must invoke this method in his ConstructProcess() 
   //  implementation in order to insures particle transportation.
   void AddTransportation();

   //Register a process to the particle type 
   // according to the ordering parameter table
   //  'true' is returned if the process is registerd successfully
   G4bool RegisterProcess(G4VProcess*            process,
			  G4ParticleDefinition*  particle);


  public:
   void UseCoupledTransportation(G4bool vl=true);

  /////////////////////////////////////////////////////////////////
  public: // with description 
   //  "SetCuts" method sets a cut value for all particle types 
   //   in the particle table
   virtual void SetCuts(); 

  public:  // with description
   //  set/get the default cut value
   //  Calling SetDefaultCutValue causes re-calcuration of cut values
   //  and physics tables just before the next event loop
   void     SetDefaultCutValue(G4double newCutValue);
   G4double GetDefaultCutValue() const;

  /////////////////////////////////////////////////////////////////////  
  public: // with description
    // Invoke BuildPhysicsTable for all processes for all particles
    // In case of "Retrieve" flag is ON, PhysicsTable will be
    // retrieved from files
    void BuildPhysicsTable();    
  
   // do PreparePhysicsTable for specified particle type
    void PreparePhysicsTable(G4ParticleDefinition* );    

   // do BuildPhysicsTable for specified particle type
    void BuildPhysicsTable(G4ParticleDefinition* );    

     // Store PhysicsTable together with both material and cut value 
    // information in files under the specified directory.
    //  (return true if files are sucessfully created)
    G4bool  StorePhysicsTable(const G4String& directory = ".");
 
    // Return true if "Retrieve" flag is ON. 
    // (i.e. PhysicsTable will be retrieved from files)
    G4bool  IsPhysicsTableRetrieved() const;
    G4bool  IsStoredInAscii() const;

    // Get directory path for physics table files.
    const G4String& GetPhysicsTableDirectory() const;

    // Set "Retrieve" flag
    // Directory path can be set together.
    // Null string (default) means directory is not changed 
    // from the current value 
    void    SetPhysicsTableRetrieved(const G4String& directory = "");
    void    SetStoredInAscii();
  
    // Reset "Retrieve" flag
    void    ResetPhysicsTableRetrieved();
    void    ResetStoredInAscii();

 ///////////////////////////////////////////////////////////////////////
  public: // with description
    // Print out the List of registered particles types
    void DumpList() const;

  public: // with description
    // Request to print out information of cut values
    // Printing will be performed when all tables are made
    void DumpCutValuesTable(G4int flag =1);

    // The following method actually trigger the print-out requested
    // by the above method. This method must be invoked by RunManager
    // at the proper moment.
    void DumpCutValuesTableIfRequested();

  public: // with description
    void  SetVerboseLevel(G4int value);
    G4int GetVerboseLevel() const;
    // set/get controle flag for output message
    //  0: Silent
    //  1: Warning message
    //  2: More

  ///////////////////////////////////////////////////////////////////////////
  public: // with description
   //  "SetCutsWithDefault" method invokes default SetCuts method
   //   Note: Cut values will not be overwriten with this method 
   //   Using default SetCuts method is recommended
   //  (i.e You do not need to implement SetCuts method) 
   void SetCutsWithDefault();   

   // Following are utility methods for SetCuts
  
   // SetCutValue sets a cut value for a particle type for the default region
   void SetCutValue(G4double aCut, const G4String& pname); 

   // GetCutValue sets a cut value for a particle type for the default region
   G4double GetCutValue(const G4String& pname) const; 

   // SetCutValue sets a cut value for a particle type for a region
   void SetCutValue(G4double aCut, const G4String& pname, const G4String& rname); 

   // Invoke SetCuts for specified particle for a region
   // If the pointer to the region is NULL, the default region is used
   // In case of "Retrieve" flag is ON, 
   // Cut values will be retrieved from files
   void SetParticleCuts(G4double cut,G4ParticleDefinition* particle,G4Region* region=0);
  void SetParticleCuts( G4double cut, const G4String& particleName, G4Region* region=0);

   // Invoke SetCuts for all particles in a region
   void SetCutsForRegion(G4double aCut, const G4String& rname);

   // Following are utility methods are obsolete
   void ResetCuts();

///////////////////////////////////////////////////////////////////
  public:   
   // Get/SetApplyCuts gets/sets the flag for ApplyCuts
   void SetApplyCuts(G4bool value, const G4String& name); 
   G4bool GetApplyCuts(const G4String& name) const; 

///////////////////////////////////////////////////////////////////////////////
  protected:  
    // do BuildPhysicsTable for make the integral schema
    void BuildIntegralPhysicsTable(G4VProcess* ,G4ParticleDefinition*  );   


  protected: 
    // Retrieve PhysicsTable from files for proccess belongng the particle.
    // Normal BuildPhysics procedure of processes will be invoked, 
    // if it fails (in case of Process's RetrievePhysicsTable returns false)
    virtual void  RetrievePhysicsTable(G4ParticleDefinition* ,	
				       const G4String& directory, 
				       G4bool          ascii = false);

   /////////////////////////////////////////////////////////////////
  protected: 
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
			   G4ProcessManager*    newManager = 0 );
 
   /////////////////////////////////////////////////////////////////
  public:
    // check consistencies of list of particles 

    void CheckParticleList();

    void DisableCheckParticleList();
 
  ////////////////////////////////////////////////////////////////////////
  protected:
   // the particle table has the complete List of existing particle types
   G4ParticleTable* theParticleTable;
   //G4ParticleTable::G4PTblDicIterator* theParticleIterator; //AND

  protected: 
   // pointer to G4UserPhysicsListMessenger
   //G4UserPhysicsListMessenger* theMessenger;

  protected:
   G4int verboseLevel;

  protected:
   // this is the default cut value for all particles
   G4double defaultCutValue;
  G4bool   isSetDefaultCutValue;

  protected:
   // pointer to ProductionCutsTable
   G4ProductionCutsTable* fCutsTable;

   // flag to determine physics table will be build from file or not
   G4bool fRetrievePhysicsTable;  
   G4bool fStoredInAscii;
 
   G4bool fIsCheckedForRetrievePhysicsTable;
   G4bool fIsRestoredCutValues;

   // directory name for physics table files 
   G4String directoryPhysicsTable;   

   // flag for displaying the range cuts & energy thresholds
   //G4int fDisplayThreshold;

  // flag for Physics Table has been built 
   //G4bool fIsPhysicsTableBuilt;

  // flag for CheckParticleList 
  G4bool fDisableCheckParticleList; 

  // PhysicsListHelper
  //G4PhysicsListHelper* thePLHelper;

  private:
   enum { FixedStringLengthForStore = 32 }; 

  //Changes for MT
  protected:
    G4int g4vuplInstanceID;
    G4RUN_DLL static G4VUPLManager subInstanceManager;
    G4ParticleTable::G4PTblDicIterator* GetParticleIterator() const;
  public:
    inline G4int GetInstanceID() const;
    static const G4VUPLManager& GetSubInstanceManager();
    //Used by Worker threads on the shared instance of
    // PL to initialize workers. Derived class re-implementing this method
    // must also call this base class method
    virtual void InitializeWorker();
    //Destroy thread-local data. Note that derived classes
    //implementing this method should still call this base class one
    virtual void TerminateWorker();
};

inline void G4VUserPhysicsList::Construct()
{
#ifdef G4VERBOSE  
  if (verboseLevel >1) G4cout << "G4VUserPhysicsList::Construct()" << G4endl;  
#endif

  InitializeProcessManager();

#ifdef G4VERBOSE  
  if (verboseLevel >1) G4cout << "Construct processes " << G4endl;  
#endif
  ConstructProcess();

}

inline G4double G4VUserPhysicsList::GetDefaultCutValue() const
{
  return defaultCutValue;
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
 G4bool  G4VUserPhysicsList::IsStoredInAscii() const
{
  return fStoredInAscii;
}

inline 
  const G4String& G4VUserPhysicsList::GetPhysicsTableDirectory() const
{
  return directoryPhysicsTable;  
}

inline 
 void  G4VUserPhysicsList::SetStoredInAscii()
{
  fStoredInAscii = true;
}
    
    
inline 
 void  G4VUserPhysicsList::ResetPhysicsTableRetrieved()
{
  fRetrievePhysicsTable = false;
  fIsRestoredCutValues = false;
  fIsCheckedForRetrievePhysicsTable=false;
}


inline 
 void  G4VUserPhysicsList::ResetStoredInAscii()
{
  fStoredInAscii = false;
}

inline
 void G4VUserPhysicsList::DisableCheckParticleList()
{
  fDisableCheckParticleList = true;
}

inline
G4int G4VUserPhysicsList::GetInstanceID() const
{
    return g4vuplInstanceID;
}

inline
const G4VUPLManager& G4VUserPhysicsList::GetSubInstanceManager()
{
    return subInstanceManager;
}
#endif

