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
// $Id: G4VUserPhysicsList.hh,v 1.41 2009-08-09 14:31:46 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// ------------------------------------------------------------
#ifndef G4VUserPhysicsList_h
#define G4VUserPhysicsList_h 1
#include "globals.hh"
#include "G4ios.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh" 
#include "G4ProductionCutsTable.hh"

class G4UserPhysicsListMessenger;
class G4VProcess;

class G4VUserPhysicsList
{
  public: 
    G4VUserPhysicsList();
    virtual ~G4VUserPhysicsList();

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
   //  !! Caution: this class must not be overriden !!
   void AddTransportation();

  /////////////////////////////////////////////////////////////////
  public: // with description 
   //  "SetCuts" method sets a cut value for all particle types 
   //   in the particle table
   virtual void SetCuts() = 0; 

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
    void DumpCutValuesTable(G4int nParticles=4);

    // The following method actually trigger the print-out requested
    // by the above method. This method must be invoked by RunManager
    // at the proper moment.
    void DumpCutValuesTableIfRequested();

    void DumpCutValues(const G4String &particle_name = "ALL");
    void DumpCutValues(G4ParticleDefinition* );

  public: // with description
    void  SetVerboseLevel(G4int value);
    G4int GetVerboseLevel() const;
    // set/get controle flag for output message
    //  0: Silent
    //  1: Warning message
    //  2: More

  ///////////////////////////////////////////////////////////////////////////
  public: // with description
   //  "SetCutsWithDefault" method sets the default cut value
   //   for all particles for the default region.
   void SetCutsWithDefault();   

   // Following are utility methods for SetCuts
  
   // SetCutValue sets a cut value for a particle type for the default region
   void SetCutValue(G4double aCut, const G4String& pname); 

   // SetCutValue sets a cut value for a particle type for a region
   void SetCutValue(G4double aCut, const G4String& pname, const G4String& rname); 

   // Invoke SetCuts for specified particle for a region
   // If the pointer to the region is NULL, the default region is used
   // In case of "Retrieve" flag is ON, 
   // Cut values will be retrieved from files
   void SetParticleCuts(G4double cut,G4ParticleDefinition* particle,G4Region* region=0);

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
 
 protected: 
 
    bool fDisableCheckParticleList;

  ////////////////////////////////////////////////////////////////////////
  protected:
    // the particle table has the complete List of existing particle types
    G4ParticleTable* theParticleTable;
    G4ParticleTable::G4PTblDicIterator* theParticleIterator;

  protected: 
  // pointer to G4UserPhysicsListMessenger
    G4UserPhysicsListMessenger* theMessenger;

  protected:
    G4int verboseLevel;

  protected:
    // this is the default cut value for all particles
    G4double defaultCutValue;

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
   G4int fDisplayThreshold;

  // flag for Physics Table has been built 
   G4bool fIsPhysicsTableBuilt;

  private:
   enum { FixedStringLengthForStore = 32 }; 


  private:
   G4bool useCoupledTransportation;
    
  public:
   inline void UseCoupledTransportation(G4bool vl=true)
   { useCoupledTransportation = vl; }

////////////////////////////////////////////////////////////////////////////
// Following method is for backward compatibility and removed soon
////////////////////////////////////////////////////////////////////////////
  protected:
  void SetCutValueForOthers(G4double) const;

};

inline
 void G4VUserPhysicsList::SetCutValueForOthers(G4double) const
{
  G4cerr << "WARNING !" << G4endl;
  G4cerr << " SetCutValueForOthers became obsolete." << G4endl;
  G4cerr << " It is harmless to remove this invokation without any side effects." << G4endl;
  G4cerr << " This dummy method implementation will be removed soon." << G4endl;
}

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
 void  G4VUserPhysicsList::DisableCheckParticleList()
{   
  fDisableCheckParticleList = true;
}


#endif

