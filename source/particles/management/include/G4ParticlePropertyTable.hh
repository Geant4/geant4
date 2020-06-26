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
// G4ParticlePropertyTable
//
// Class description:
//
// This class manages properties of a particle which are properties
// in G4ParticlePropertyTable class. This class is a singleton.

// Author: H.Kurashige, 9 June 2003 - First implementation
// --------------------------------------------------------------------
#ifndef G4ParticlePropertyTable_hh
#define G4ParticlePropertyTable_hh 1

#include <vector>

#include "globals.hh"
#include "G4ios.hh"

#include "G4ParticlePropertyData.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
 
class G4ParticlePropertyTable
{
  public:

    G4ParticlePropertyTable(const G4ParticlePropertyTable&) = delete;  
    G4ParticlePropertyTable& operator=(const G4ParticlePropertyTable&) = delete;

   ~G4ParticlePropertyTable();

    static G4ParticlePropertyTable* GetParticlePropertyTable();
      // Return the pointer to G4ParticlePropertyTable object
      // G4ParticlePropertyTable is a "singleton" and can get its pointer 
      // by this function. At the first time of calling this function, 
      // the G4ParticleTable object is instantiated 

    G4ParticlePropertyData* GetParticleProperty(const G4String& aParticleName);
    G4ParticlePropertyData* GetParticleProperty(const G4ParticleDefinition* aP);
      // Return the pointer to G4ParticlePropertyData object,
      // which contains properties for the particle specified.
      // (return 0 if the specified particle does not exist)
 
    G4bool SetParticleProperty(const G4ParticlePropertyData& newProperty);
      // Change particle properties for the particle specified.
      // Return true if properties are successfully set 

    void Clear();
      // Clear and destroy data

    void  SetVerboseLevel(G4int value);
    G4int GetVerboseLevel() const;
      // Control flag for output message
      //  0: Silent
      //  1: Warning message
      //  2: More

  protected: 

    G4ParticlePropertyTable();
      // Hidden default constructor; this class is a singleton  
  
    G4ParticleTable* fParticleTable = nullptr;  

    std::vector<G4ParticlePropertyData*> arrayDataObject; 
  
  private:

    G4int verboseLevel = 1;
    static G4ThreadLocal G4ParticlePropertyTable* fgParticlePropertyTable;
};

#endif
