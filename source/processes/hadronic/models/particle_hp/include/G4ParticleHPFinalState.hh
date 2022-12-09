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
// 080721 Create adjust_final_state method by T. Koi  
// 080801 Introduce theNDLDataA,Z which has A and Z of NDL data by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPFinalState_h
#define G4ParticleHPFinalState_h

#include "G4ParticleHPManager.hh"
#include "G4Material.hh"
#include "G4HadFinalState.hh"
#include "G4ParticleHPNames.hh"
#include "G4ParticleHPVector.hh"
#include "G4HadProjectile.hh"
#include "G4Neutron.hh"
#include "G4Cache.hh"

class G4ParticleDefinition;

class G4ParticleHPFinalState
{
public:

  G4ParticleHPFinalState()
  { 
    hasFSData = true; 
    hasXsec = true;
    hasAnyData = true;
    theBaseZ = 0;
    theBaseA = 0;
    theBaseM = 0;

    theNDLDataZ = 0;
    theNDLDataA = 0;
    theNDLDataM = 0;

    theProjectile = G4Neutron::Neutron();

    theResult.Put( 0 );

    secID = -1;
  }
  
  virtual ~G4ParticleHPFinalState()
  {
    if (theResult.Get() != 0) delete theResult.Get();
  }

  void Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType,
             G4ParticleDefinition* projectile)
  {
    G4int M = 0;
    Init ( A, Z, M, dirName, aFSType,const_cast<G4ParticleDefinition*>(projectile));
  }
  virtual void Init (G4double A, G4double Z, G4int M, G4String & dirName,
                     G4String & aFSType, G4ParticleDefinition* ) = 0;
  virtual G4HadFinalState * ApplyYourself(const G4HadProjectile & ) 
  {
    throw G4HadronicException(__FILE__, __LINE__, "G4HadFinalState * ApplyYourself(const G4HadProjectile & theTrack) needs implementation");
    return 0;
  }
  
  // this would better be done templating G4ParticleHPChannel..., 
  // 
  virtual G4ParticleHPFinalState * New() = 0;
  
  G4bool HasXsec() { return hasXsec; }
  G4bool HasFSData() { return hasFSData; }
  G4bool HasAnyData() { return hasAnyData; }
  
  virtual G4double GetXsec(G4double ) { return 0; }
  virtual G4ParticleHPVector * GetXsec() { return 0; }
  
  void     SetA_Z(G4double anA, G4double aZ, G4int aM=0) {theBaseA = anA; theBaseZ = aZ; theBaseM=aM; }
  G4double GetZ() { return theBaseZ; }
  G4double GetN() { return theBaseA; }
  G4double GetA() { return theBaseA; }
  G4int GetM() { return theBaseM; }

  void SetAZMs(G4double anA, G4double aZ, G4int aM, G4ParticleHPDataUsed used) 
  {
    theBaseA = anA; theBaseZ = aZ; theBaseM=aM; 
    theNDLDataA=(G4int)used.GetA();
    theNDLDataZ=(G4int)used.GetZ();
    theNDLDataM=used.GetM();
  }
  
  void SetProjectile( G4ParticleDefinition* projectile )
  {
    theProjectile = projectile;
  }

protected:
  
    void adjust_final_state ( G4LorentzVector );
  
    G4bool hasXsec;
    G4bool hasFSData;
    G4bool hasAnyData;
    G4ParticleHPNames theNames;
  
    G4Cache< G4HadFinalState* > theResult;
    G4ParticleDefinition* theProjectile;
  
    G4double theBaseA;
    G4double theBaseZ;
    G4int theBaseM;

    G4int theNDLDataZ;
    G4int theNDLDataA;
    G4int theNDLDataM;

    G4int secID;  // Creator model ID for the secondaries created by this class or derived ones
};
#endif
