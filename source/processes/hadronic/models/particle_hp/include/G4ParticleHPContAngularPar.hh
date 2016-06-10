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
//
// 080718 Add ClearHistories method and related class member
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPContAngularPar_h
#define G4ParticleHPContAngularPar_h 1

#include "G4ios.hh"
#include <fstream>
#include "globals.hh"
#include "G4ParticleHPList.hh"
#include "G4ReactionProduct.hh"
#include "G4ParticleHPInterpolator.hh"
#include "G4InterpolationManager.hh"
#include <set>
class G4ParticleDefinition;

class G4ParticleHPContAngularPar
{
  public:
  
  G4ParticleHPContAngularPar()
  {
    theAngular = 0;
    currentMeanEnergy = -2;
    fresh = true;
    theMinEner = DBL_MAX;
    theMaxEner = -DBL_MAX;
  }

  G4ParticleHPContAngularPar(G4ParticleDefinition* projectile); 

  ~G4ParticleHPContAngularPar()
  {
    if(theAngular!=0) delete [] theAngular;
  }
  
  void Init(std::istream & aDataFile, G4ParticleDefinition* projectile);
  
  G4ReactionProduct * Sample(G4double anEnergy, G4double massCode, G4double mass, 
                             G4int angularRep, G4int interpol);
  
  G4double GetEnergy() { 
    if( getenv("G4PHPTEST") ) G4cout << this << " G4ParticleHPContAngularPar::GetEnergy " << theEnergy <<  " nE " << nEnergies << G4endl;
    return theEnergy; }
  
  void SetPrimary(G4ReactionProduct * aPrimary)
  {
    thePrimary = aPrimary;
  }
  
  void SetTarget(G4ReactionProduct * aTarget)
  {
    theTarget = aTarget;
  }
  
  void SetTargetCode(G4double aTargetCode) { theTargetCode = aTargetCode; }
  
  void SetInterpolation(G4int theInterpolation)
  {
    theManager.Init(theInterpolation, nEnergies); // one range only
  }

  void BuildByInterpolation(G4double anEnergy, G4InterpolationScheme aScheme, 
             G4ParticleHPContAngularPar & store1, 
             G4ParticleHPContAngularPar & store2); // hmmmm, this interpolates legendre coefficients. Dangerous @@@

  void PrepareTableInterpolation(const G4ParticleHPContAngularPar* angularPrev);
  
  G4double MeanEnergyOfThisInteraction()
  {
    G4double result;
    if(currentMeanEnergy<-1)
    {
      return 0;
      // throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPContAngularPar: Logical error in Product class");
    }
    else
    {
      result = currentMeanEnergy;
    }
    currentMeanEnergy = -2;
    return result;
  }
  
  G4int GetNEnergies() const 
  {
    return nEnergies; 
  }
  G4int GetNDiscreteEnergies() const 
  {
    return nDiscreteEnergies; 
  }
  std::set<G4double> GetEnergiesTransformed() const 
  {
    return theEnergiesTransformed;
  }
  G4int GetNEnergiesTransformed() const 
  {
    return theEnergiesTransformed.size();
  }
  G4double GetMinEner() const 
  {
    return theMinEner;
  }
  G4double GetMaxEner() const 
  { 
    return theMaxEner;
  }
  std::map<G4double,G4int> GetDiscreteEnergiesOwn() const 
  {
    return theDiscreteEnergiesOwn;
  }
  G4ParticleHPList* GetAngDataList() const {
    return theAngular; 
  }
  
private:
  
  // incoming particle
  G4double theEnergy; 
  
  // number of exit channel energies
  G4int nEnergies; 
  // number of discrete exit channels
  G4int nDiscreteEnergies;
  // number of angular paramerers per channel
  G4int nAngularParameters;
  // knows the interpolation between List labels
  G4InterpolationManager theManager; 
  // on per exit-channel energy
  G4ParticleHPList * theAngular; 
  
  G4ParticleHPInterpolator theInt;
  
  G4double theTargetCode;
  G4ReactionProduct * theTarget;
  G4ReactionProduct * thePrimary;
  
  G4double currentMeanEnergy;

//080718
   public:
      void ClearHistories(){ fresh = true; };

  void Dump();
   private:
      G4bool fresh; 
      G4double remaining_energy; // represent energy rest of cascade chain

   G4ParticleDefinition* theProjectile;

  G4bool adjustResult; // if not set it will not force the conservation of energy in angularRep==1, but will sample the particle energy according to the database

  G4double theMinEner;
  G4double theMaxEner;
  std::set<G4double> theEnergiesTransformed;
  std::set<G4double> theDiscreteEnergies;
  std::map<G4double,G4int> theDiscreteEnergiesOwn;

};
#endif
