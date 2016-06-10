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
// $Id: G4ParticlePropertyTable.hh 67971 2013-03-13 10:13:24Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: 
// ---------------- G4ParticlePropertyTable ----------------
// first implementation by H Kurashige 9 June 2003
// ------------------------------------------------------------

#ifndef G4ParticlePropertyTable_h
#define G4ParticlePropertyTable_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "G4ParticlePropertyData.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
 
class G4ParticlePropertyTable
{
  // Class Description
  //  This class manages properties of a particle which are 
  //  properties in G4ParticlePropertyTable class.
  //  This class is a singleton.
 
 protected: 
  // hide default constructor because this class is a singleton  
  G4ParticlePropertyTable();
  G4ParticlePropertyTable(const G4ParticlePropertyTable &right);
      
  G4ParticlePropertyTable & operator=(const G4ParticlePropertyTable &right);
  
 public:
  G4int operator==(const G4ParticlePropertyTable &right) const;
  G4int operator!=(const G4ParticlePropertyTable &right) const;

 public:
  virtual ~G4ParticlePropertyTable();

 public: //With Description
  static G4ParticlePropertyTable* GetParticlePropertyTable();
   // return the pointer to G4ParticlePropertyTable object
   // G4ParticlePropertyTable is a "singleton" and can get its pointer 
   // by this function. At the first time of calling this function, 
   // the G4ParticleTable object is instantiated 

  G4ParticlePropertyData* GetParticleProperty(const G4String& aParticleName);
  G4ParticlePropertyData* GetParticleProperty(const G4ParticleDefinition* aParticle);
  // return the pointer to  G4ParticlePropertyData object,
  // which contains properties for the particle specified.
  // (return 0 if the specified particle does not exist)
 
  G4bool SetParticleProperty(const G4ParticlePropertyData& newProperty);
  // change particle properties for the particle specified.
  // return true if properties are sucessfully set 

  void   Clear();
  // clear and destroy arrayDataObject

 public:
  void  SetVerboseLevel(G4int value);
  G4int GetVerboseLevel() const;
  // controle flag for output message
  //  0: Silent
  //  1: Warning message
  //  2: More

 protected:
  G4ParticleTable* fParticleTable;  
  
 private:
  G4int verboseLevel;
  static G4ThreadLocal G4ParticlePropertyTable*  fgParticlePropertyTable;
  
 protected:
  std::vector<G4ParticlePropertyData*> arrayDataObject; 
};


#endif
