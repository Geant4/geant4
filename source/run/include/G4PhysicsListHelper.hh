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
// G4PhysicsListHelper
//
// Class description:
//
// Helper class for physics lists, to register processes according
// to the ordering parameter table. This class is a singleton.

// Author: H.Kurashige, 29 April 2011
// --------------------------------------------------------------------
#ifndef G4PhysicsListHelper_hh
#define G4PhysicsListHelper_hh 1

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicsListOrderingParameter.hh"
#include "G4ThreadLocalSingleton.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <vector>

class G4VProcess;

class G4PhysicsListHelper
{
    friend class G4ThreadLocalSingleton<G4PhysicsListHelper>;

  public:
    // Returns the pointer to the physics list helper
    static G4PhysicsListHelper* GetPhysicsListHelper();

    // Registers a process to the particle type according to the ordering
    // parameter table. Returns 'true' if process is successfully registered.
    G4bool RegisterProcess(G4VProcess* process, G4ParticleDefinition* particle);

    // User must invoke this method in his ConstructProcess() implementation
    // in order to enable particle transportation.
    void AddTransportation();

    // Set flag for using G4CoupledTransportation.
    void UseCoupledTransportation(G4bool vl = true);

    // Change the thresholds for killing looping tracks in transportation.
    void UseHighLooperThresholds() { theLooperThresholds = 2; }
    void UseLowLooperThresholds() { theLooperThresholds = 0; }

    // Check consistencies of list of particles.
    void CheckParticleList() const;

    // Dump OrdingParameterTable.
    void DumpOrdingParameterTable(G4int subType = -1) const;

    G4PhysicsListOrderingParameter GetOrdingParameter(G4int subType) const;

    // set/get controle flag for output message
    //  0: Silent
    //  1: Warning message
    //  2: More
    void SetVerboseLevel(G4int value);
    G4int GetVerboseLevel() const;

  private:
    // Hidden constructor and destructor.
    G4PhysicsListHelper();
    ~G4PhysicsListHelper();

    void ReadOrdingParameterTable();
    void ReadInDefaultOrderingParameter();

  private:
    using G4OrdParamTable = std::vector<G4PhysicsListOrderingParameter>;

    static G4ThreadLocal G4PhysicsListHelper* pPLHelper;

    // The particle table has the complete List of existing particle types.
    G4ParticleTable* theParticleTable = nullptr;
    G4ParticleTable::G4PTblDicIterator* aParticleIterator = nullptr;

    G4bool useCoupledTransportation = false;
    G4int theLooperThresholds = 1;  //  0 = Low,  1 = default, 2 = high
    G4VProcess* theTransportationProcess = nullptr;

    G4int verboseLevel = 1;

    G4OrdParamTable* theTable = nullptr;
    G4int sizeOfTable = 0;
    G4String ordParamFileName = "";
};

// Inline methods implementations

inline void G4PhysicsListHelper::UseCoupledTransportation(G4bool vl)
{
  useCoupledTransportation = vl;
}

inline void G4PhysicsListHelper::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline G4int G4PhysicsListHelper::GetVerboseLevel() const
{
  return verboseLevel;
}

#endif
