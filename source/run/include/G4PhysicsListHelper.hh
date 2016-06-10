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
// ------------------------------------------------------------
//	GEANT 4 class header file 
// Class Description:
//      This class is a helper class for physics lists to register processes 
//      according to the ordering parameter table 
//      This class is a singleton
// ------------------------------------------- 
//	History
//        first version                   29 Apr. 2011 by H.Kurashige 
// ------------------------------------------------------------

#ifndef G4PhysicsListHelper_h
#define G4PhysicsListHelper_h 1
#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh" 
#include "G4PhysicsListOrderingParameter.hh" 

class G4VProcess;

class G4PhysicsListHelper
{
  private:
   // Hide constructor and destructor 
   G4PhysicsListHelper();
   virtual ~G4PhysicsListHelper();

  public:  // with description
   // This method gives the ponter to the physics list helper 
   static G4PhysicsListHelper* GetPhysicsListHelper(); 
  
   //Register a process to the particle type 
   // according to the ordering parameter table
   //  'true' is returned if the process is registerd successfully
   G4bool RegisterProcess(G4VProcess*            process,
			  G4ParticleDefinition*  particle);

   //  User must invoke this method in his ConstructProcess() 
   //  implementation in order to insures particle transportation.
   void AddTransportation();
   
   //  Set flag for using CoupledTransportation
   void UseCoupledTransportation(G4bool vl=true);
 
  /////////////////////////////////////////////////////////////////
  public:
    // check consistencies of list of particles 
    void CheckParticleList() const;

  ///////////////////////////////////////////////////////////////////////
  public: 
  // Dump OrdingParameterTable
    void DumpOrdingParameterTable(G4int subType = -1) const;
    G4PhysicsListOrderingParameter GetOrdingParameter(G4int subType) const;

  private: 
    void ReadOrdingParameterTable();
    void ReadInDefaultOrderingParameter();

  ///////////////////////////////////////////////////////////////////////
  public: // with description
    void  SetVerboseLevel(G4int value);
    G4int GetVerboseLevel() const;
    // set/get controle flag for output message
    //  0: Silent
    //  1: Warning message
    //  2: More

  ////////////////////////////////////////////////////////////////////////
  private:
    static G4ThreadLocal G4PhysicsListHelper* pPLHelper;

    // the particle table has the complete List of existing particle types
    G4ParticleTable* theParticleTable;
    G4ParticleTable::G4PTblDicIterator* aParticleIterator;

    G4bool useCoupledTransportation;
    G4VProcess* theTransportationProcess;
 
    G4int verboseLevel;

  private:
    typedef std::vector<G4PhysicsListOrderingParameter> G4OrdParamTable;
    G4OrdParamTable* theTable;
    G4int            sizeOfTable;
    G4String         ordParamFileName;
};

 
inline 
void G4PhysicsListHelper::UseCoupledTransportation(G4bool vl)
{ 
  useCoupledTransportation = vl; 
}

inline
 void  G4PhysicsListHelper::SetVerboseLevel(G4int value)
{  
  verboseLevel = value;
}
    
inline
  G4int G4PhysicsListHelper::GetVerboseLevel() const
{
  return verboseLevel;
}


#endif

