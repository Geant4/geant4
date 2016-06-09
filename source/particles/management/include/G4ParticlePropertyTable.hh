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
// $Id: G4ParticlePropertyTable.hh,v 1.3 2003/12/09 15:35:46 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
      
  const G4ParticlePropertyTable & operator=(const G4ParticlePropertyTable &right);
  
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
  static G4ParticlePropertyTable*  fgParticlePropertyTable;
  
 protected:
  std::vector<G4ParticlePropertyData*> arrayDataObject; 
};


#endif









